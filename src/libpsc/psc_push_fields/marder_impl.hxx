
#pragma once

#include "marder.hxx"
#include "fields.hxx"
#include "writer_mrc.hxx"
#include "../libpsc/psc_output_fields/fields_item_fields.hxx"
#include "../libpsc/psc_output_fields/fields_item_moments_1st.hxx"
#include "../libpsc/psc_bnd/psc_bnd_impl.hxx"

#include <gtensor/reductions.h>
#include <mrc_io.h>

#ifdef USE_CUDA
void cuda_marder_correct_yz(MfieldsStateCuda& mflds, MfieldsCuda& mf, int p,
                            Float3 fac, Int3 ly, Int3 ry, Int3 lz, Int3 rz);
void cuda_marder_correct_xyz(MfieldsStateCuda& mflds, MfieldsCuda& mf, int p,
                             Float3 fac, Int3 lx, Int3 rx, Int3 ly, Int3 ry,
                             Int3 lz, Int3 rz);
#endif

namespace psc
{
namespace marder
{

// ----------------------------------------------------------------------
// correct
//
// Do the modified marder correction (See eq.(5, 7, 9, 10) in Mardahl and
// Verboncoeur, CPC, 1997)

#define define_dxdydz(dx, dy, dz)                                              \
  int dx _mrc_unused = (grid.isInvar(0)) ? 0 : 1;                              \
  int dy _mrc_unused = (grid.isInvar(1)) ? 0 : 1;                              \
  int dz _mrc_unused = (grid.isInvar(2)) ? 0 : 1

#define psc_foreach_3d_more(psc, p, ix, iy, iz, l, r)                          \
  {                                                                            \
    int __ilo[3] = {-l[0], -l[1], -l[2]};                                      \
    int __ihi[3] = {grid.ldims[0] + r[0], grid.ldims[1] + r[1],                \
                    grid.ldims[2] + r[2]};                                     \
    for (int iz = __ilo[2]; iz < __ihi[2]; iz++) {                             \
      for (int iy = __ilo[1]; iy < __ihi[1]; iy++) {                           \
        for (int ix = __ilo[0]; ix < __ihi[0]; ix++)

#define psc_foreach_3d_more_end                                                \
  }                                                                            \
  }                                                                            \
  }

template <typename MfieldsState, typename Mfields>
inline void correct(MfieldsState& mflds, Mfields& mf,
                    typename MfieldsState::real_t diffusion)
{
  const auto& grid = mflds.grid();
  define_dxdydz(dx, dy, dz);

  // FIXME: how to choose diffusion parameter properly?
  double deltax = grid.domain.dx[0]; // FIXME double/float
  double deltay = grid.domain.dx[1];
  double deltaz = grid.domain.dx[2];

  for (int p = 0; p < mf.n_patches(); p++) {
    int l_cc[3] = {0, 0, 0}, r_cc[3] = {0, 0, 0};
    int l_nc[3] = {0, 0, 0}, r_nc[3] = {0, 0, 0};
    for (int d = 0; d < 3; d++) {
      if (grid.bc.fld_lo[d] == BND_FLD_CONDUCTING_WALL &&
          grid.atBoundaryLo(p, d)) {
        l_cc[d] = -1;
        l_nc[d] = -1;
      }
      if (grid.bc.fld_hi[d] == BND_FLD_CONDUCTING_WALL &&
          grid.atBoundaryHi(p, d)) {
        r_cc[d] = -1;
        r_nc[d] = 0;
      }
    }

    auto flds_ = make_Fields3d<dim_xyz>(mflds[p]);
    auto f_ = make_Fields3d<dim_xyz>(mf[p]);
    if (!grid.isInvar(0)) {
      int l[3] = {l_cc[0], l_nc[1], l_nc[2]};
      int r[3] = {r_cc[0], r_nc[1], r_nc[2]};
      psc_foreach_3d_more(ppsc, p, ix, iy, iz, l, r)
      {
        flds_(EX, ix, iy, iz) += (f_(0, ix + dx, iy, iz) - f_(0, ix, iy, iz)) *
                                 .5 * grid.dt * diffusion / deltax;
      }
      psc_foreach_3d_more_end;
    }

    {
      int l[3] = {l_nc[0], l_cc[1], l_nc[2]};
      int r[3] = {r_nc[0], r_cc[1], r_nc[2]};
      psc_foreach_3d_more(ppsc, p, ix, iy, iz, l, r)
      {
        flds_(EY, ix, iy, iz) += (f_(0, ix, iy + dy, iz) - f_(0, ix, iy, iz)) *
                                 .5 * grid.dt * diffusion / deltay;
      }
      psc_foreach_3d_more_end;
    }

    {
      int l[3] = {l_nc[0], l_nc[1], l_cc[2]};
      int r[3] = {r_nc[0], r_nc[1], r_cc[2]};
      psc_foreach_3d_more(ppsc, p, ix, iy, iz, l, r)
      {
        flds_(EZ, ix, iy, iz) += (f_(0, ix, iy, iz + dz) - f_(0, ix, iy, iz)) *
                                 .5 * grid.dt * diffusion / deltaz;
      }
      psc_foreach_3d_more_end;
    }
  }
}

#undef psc_foreach_3d_more
#undef psc_foreach_3d_more_end

#ifdef USE_CUDA

inline void correct(MfieldsStateCuda& mflds, MfieldsCuda& mf, float diffusion)
{
  const auto& grid = mflds.grid();

  Float3 fac;
  fac[0] = .5 * grid.dt * diffusion / grid.domain.dx[0];
  fac[1] = .5 * grid.dt * diffusion / grid.domain.dx[1];
  fac[2] = .5 * grid.dt * diffusion / grid.domain.dx[2];

  // OPT, do all patches in one kernel
  for (int p = 0; p < mflds.n_patches(); p++) {
    int l_cc[3] = {0, 0, 0}, r_cc[3] = {0, 0, 0};
    int l_nc[3] = {0, 0, 0}, r_nc[3] = {0, 0, 0};
    for (int d = 0; d < 3; d++) {
      if (grid.bc.fld_lo[d] == BND_FLD_CONDUCTING_WALL &&
          grid.atBoundaryLo(p, d)) {
        l_cc[d] = -1;
        l_nc[d] = -1;
      }
      if (grid.bc.fld_hi[d] == BND_FLD_CONDUCTING_WALL &&
          grid.atBoundaryHi(p, d)) {
        r_cc[d] = -1;
        r_nc[d] = 0;
      }
    }

    Int3 ldims = grid.ldims;

    Int3 lx = {l_cc[0], l_nc[1], l_nc[2]};
    Int3 rx = {r_cc[0] + ldims[0], r_nc[1] + ldims[1], r_nc[2] + ldims[2]};

    Int3 ly = {l_nc[0], l_cc[1], l_nc[2]};
    Int3 ry = {r_nc[0] + ldims[0], r_cc[1] + ldims[1], r_nc[2] + ldims[2]};

    Int3 lz = {l_nc[0], l_nc[1], l_cc[2]};
    Int3 rz = {r_nc[0] + ldims[0], r_nc[1] + ldims[1], r_cc[2] + ldims[2]};

    if (grid.isInvar(0)) {
      cuda_marder_correct_yz(mflds, mf, p, fac, ly, ry, lz, rz);
    } else {
      cuda_marder_correct_xyz(mflds, mf, p, fac, lx, rx, ly, ry, lz, rz);
    }
  }
}
#endif

} // namespace marder
} // namespace psc

template <typename MP, typename MFS, typename MF, typename D, typename ITEM_RHO,
          typename BND>
class MarderCommon
{
public:
  using Mparticles = MP;
  using MfieldsState = MFS;
  using Mfields = MF;
  using dim_t = D;
  using Item_rho_t = ITEM_RHO;
  using Bnd = BND;
  using real_t = typename Mfields::real_t;

  // FIXME: checkpointing won't properly restore state
  // FIXME: if the subclass creates objects, it'd be cleaner to have them
  // be part of the subclass

  MarderCommon(const Grid_t& grid, real_t diffusion, int loop, bool dump)
    : grid_{grid},
      diffusion_{diffusion},
      loop_{loop},
      dump_{dump},
      bnd_{grid, grid.ibn},
      rho_{grid, 1, grid.ibn},
      res_{grid, 1, grid.ibn}
  {
    if (dump_) {
      io_.open("marder");
    }
  }

  // ----------------------------------------------------------------------
  // print_max

  static void print_max(Mfields& mf)
  {
    real_t max_err = gt::norm_linf(mf.storage());
    MPI_Allreduce(MPI_IN_PLACE, &max_err, 1,
                  Mfields_traits<Mfields>::mpi_dtype(), MPI_MAX,
                  mf.grid().comm());
    mpi_printf(mf.grid().comm(), "marder: err %g\n", max_err);
  }

  // ----------------------------------------------------------------------
  // calc_aid_fields

  template <typename E>
  void calc_aid_fields(MfieldsState& mflds, const E& rho)
  {
    const auto& grid = mflds.grid();
    auto item_dive = Item_dive<MfieldsState>{};
    auto dive = psc::mflds::interior(grid, item_dive(mflds));

    if (dump_) {
      static int cnt;
      io_.begin_step(cnt, cnt); // ppsc->timestep, ppsc->timestep * ppsc->dt);
      cnt++;
      io_.write(rho, grid, "rho", {"rho"});
      io_.write(dive, grid, "dive", {"dive"});
      io_.end_step();
    }

    auto& res = res_.storage();

    psc::mflds::interior(grid, res) = dive - rho;
    // FIXME, why is this necessary?
    bnd_.fill_ghosts(res_, 0, 1);
  }

  // ----------------------------------------------------------------------
  // correct
  //
  // Do the modified marder correction (See eq.(5, 7, 9, 10) in Mardahl and
  // Verboncoeur, CPC, 1997)

  void correct(MfieldsState& mflds)
  {
    auto& grid = mflds.grid();
    // FIXME: how to choose diffusion parameter properly?

    double inv_sum = 0.;
    for (int d = 0; d < 3; d++) {
      if (!grid.isInvar(d)) {
        inv_sum += 1. / sqr(grid.domain.dx[d]);
      }
    }
    double diffusion_max = 1. / 2. / (.5 * grid.dt) / inv_sum;
    double diffusion = diffusion_max * diffusion_;
    psc::marder::correct(mflds, res_, diffusion);
  }

  // private:
  const Grid_t& grid_;
  real_t diffusion_; //< diffusion coefficient for Marder correction
  int loop_;         //< execute this many relaxation steps in a loop
  bool dump_;        //< dump div_E, rho
  Bnd bnd_;
  Mfields rho_;
  Mfields res_;
  WriterMRC io_; //< for debug dumping
};

template <typename _Mparticles, typename _MfieldsState, typename _Mfields,
          typename D>
struct Marder_
  : public MarderCommon<_Mparticles, _MfieldsState, _Mfields, D,
                        Moment_rho_1st_nc<typename _Mfields::Storage, D>, Bnd_>
{
  using Base =
    MarderCommon<_Mparticles, _MfieldsState, _Mfields, D,
                 Moment_rho_1st_nc<typename _Mfields::Storage, D>, Bnd_>;
  using Mparticles = typename Base::Mparticles;
  using MfieldsState = typename Base::MfieldsState;
  using Mfields = typename Base::Mfields;
  using dim_t = typename Base::dim_t;
  using Bnd = typename Base::Bnd;
  using real_t = typename Mfields::real_t;
  using Moment_t = typename Base::Item_rho_t;
  using Base::bnd_;
  using Base::diffusion_;
  using Base::dump_;
  using Base::io_;
  using Base::loop_;
  using Base::res_;
  using Base::rho_;

  using Base::Base;

  // ----------------------------------------------------------------------
  // operator()

  void operator()(MfieldsState& mflds, Mparticles& mprts)
  {
    const auto& grid = mflds.grid();
    auto rho = Moment_t{grid};
    rho_.storage() = psc::mflds::interior(grid, rho(mprts));

    // need to fill ghost cells first (should be unnecessary with only variant
    // 1) FIXME
    bnd_.fill_ghosts(mflds, EX, EX + 3);

    for (int i = 0; i < loop_; i++) {
      Base::calc_aid_fields(mflds, rho_.storage());
      Base::print_max(res_);
      Base::correct(mflds);
      bnd_.fill_ghosts(mflds, EX, EX + 3);
    }
  }
};

#undef define_dxdydz
