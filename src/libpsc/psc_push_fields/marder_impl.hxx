
#pragma once

#include "marder.hxx"
#include "fields.hxx"
#include "writer_mrc.hxx"
#include "../libpsc/psc_output_fields/fields_item_fields.hxx"
#include "../libpsc/psc_output_fields/psc_output_fields_item_moments_1st_nc.cxx"
#include "../libpsc/psc_bnd/psc_bnd_impl.hxx"

#include <mrc_io.h>

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

} // namespace marder
} // namespace psc

template <typename _Mparticles, typename _MfieldsState, typename _Mfields,
          typename D>
struct Marder_ : MarderBase
{
  using Mparticles = _Mparticles;
  using MfieldsState = _MfieldsState;
  using Mfields = _Mfields;
  using dim_t = D;
  using real_t = typename Mfields::real_t;
  using fields_view_t = typename Mfields::fields_view_t;
  using Moment_t = Moment_rho_1st_nc<Mfields, dim_t>;

  Marder_(const Grid_t& grid, real_t diffusion, int loop, bool dump)
    : grid_{grid},
      diffusion_{diffusion},
      loop_{loop},
      dump_{dump},
      bnd_{grid, grid.ibn},
      bnd_mf_{grid, grid.ibn},
      rho_{grid, 1, grid.ibn},
      res_{grid, 1, grid.ibn}
  {
    if (dump_) {
      io_.open("marder");
    }
  }

  // FIXME: checkpointing won't properly restore state
  // FIXME: if the subclass creates objects, it'd be cleaner to have them
  // be part of the subclass

  // ----------------------------------------------------------------------
  // calc_aid_fields

  void calc_aid_fields(MfieldsState& mflds)
  {
    auto dive = Item_dive<MfieldsState>(mflds);

    if (dump_) {
      static int cnt;
      io_.begin_step(cnt, cnt); // ppsc->timestep, ppsc->timestep * ppsc->dt);
      cnt++;
      io_.write(rho_.gt(), rho_.grid(), "rho", {"rho"});
      {
        io_.write(psc::mflds::interior(dive.grid(), dive.gt()), dive.grid(),
                  "dive", {"dive"});
      }
      io_.end_step();
    }

    // res_.assign(dive);
    auto&& dive_gt = psc::mflds::interior(dive.grid(), dive.gt());
    for (int p = 0; p < res_.n_patches(); p++) {
      for (int m = 0; m < res_.n_comps(); m++) {
        kg::Box3 box = res_.box();
        Int3 ijk;
        for (ijk[2] = 0; ijk[2] < 2 * box.ib(2) + box.im(2); ijk[2]++) {
          for (ijk[1] = 0; ijk[1] < 2 * box.ib(1) + box.im(1); ijk[1]++) {
            for (ijk[0] = 0; ijk[0] < 2 * box.ib(0) + box.im(0); ijk[0]++) {
              res_(m, ijk[0], ijk[1], ijk[2], p) =
                dive_gt(ijk[0], ijk[1], ijk[2], m, p);
            }
          }
        }
      }
    }

    Int3 bnd = res_.ibn();
    res_.storage().view(_s(bnd[0], -bnd[0]), _s(bnd[1], -bnd[1]),
                        _s(bnd[2], -bnd[2])) =
      res_.storage().view(_s(bnd[0], -bnd[0]), _s(bnd[1], -bnd[1]),
                          _s(bnd[2], -bnd[2])) -
      rho_.storage().view(_all, _all, _all);
    // FIXME, why is this necessary?
    bnd_mf_.fill_ghosts(res_, 0, 1);
  }

  // ----------------------------------------------------------------------
  // max

  static void print_max(Mfields& mf)
  {
    real_t max_err = 0.;
    for (int p = 0; p < mf.n_patches(); p++) {
      auto f_ = make_Fields3d<dim_xyz>(mf[p]);
      mf.Foreach_3d(0, 0, [&](int i, int j, int k) {
        max_err = std::max(max_err, std::abs(f_(0, i, j, k)));
      });
    }
    MPI_Allreduce(MPI_IN_PLACE, &max_err, 1,
                  Mfields_traits<Mfields>::mpi_dtype(), MPI_MAX,
                  mf.grid().comm());
    mpi_printf(mf.grid().comm(), "marder: err %g\n", max_err);
  }

  // ----------------------------------------------------------------------
  // correct

  void correct(MfieldsState& mflds)
  {
    auto& grid = mflds.grid();

    double inv_sum = 0.;
    int nr_levels;
    for (int d = 0; d < 3; d++) {
      if (!grid.isInvar(d)) {
        inv_sum += 1. / sqr(grid.domain.dx[d]);
      }
    }
    double diffusion_max = 1. / 2. / (.5 * grid.dt) / inv_sum;
    double diffusion = diffusion_max * diffusion_;

    psc::marder::correct(mflds, res_, diffusion);
  }

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
      calc_aid_fields(mflds);
      print_max(res_);
      correct(mflds);
      bnd_.fill_ghosts(mflds, EX, EX + 3);
    }
  }

  // private:
  real_t diffusion_; //< diffusion coefficient for Marder correction
  int loop_;         //< execute this many relaxation steps in a loop
  bool dump_;        //< dump div_E, rho

  const Grid_t& grid_;
  Bnd_<MfieldsState> bnd_;
  Bnd_<Mfields> bnd_mf_;
  Mfields rho_;
  Mfields res_;
  WriterMRC io_; //< for debug dumping
};

#undef define_dxdydz
