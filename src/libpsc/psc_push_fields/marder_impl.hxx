
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
void cuda_marder_correct_yz(MfieldsStateCuda::Storage& gt_mflds,
                            MfieldsCuda::Storage& gt_mf, int p, Float3 fac,
                            Int3 ly, Int3 ry, Int3 lz, Int3 rz);
void cuda_marder_correct_xyz(MfieldsStateCuda::Storage& gt_mflds,
                             MfieldsCuda::Storage& gt_mf, int p, Float3 fac,
                             Int3 lx, Int3 rx, Int3 ly, Int3 ry, Int3 lz,
                             Int3 rz);
#endif

namespace psc
{
namespace marder
{

namespace detail
{

inline void find_limits(const Grid_t& grid, int p, Int3& lx, Int3& rx, Int3& ly,
                        Int3& ry, Int3& lz, Int3& rz)
{
  Int3 l_cc = {0, 0, 0}, r_cc = {0, 0, 0};
  Int3 l_nc = {0, 0, 0}, r_nc = {0, 0, 0};
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
  // FIXME, for conducting wall the signs here need checking...
  lx = -Int3{l_cc[0], l_nc[1], l_nc[2]} + grid.ibn;
  rx = Int3{r_cc[0], r_nc[1], r_nc[2]} + grid.ldims + grid.ibn;

  ly = -Int3{l_nc[0], l_cc[1], l_nc[2]} + grid.ibn;
  ry = Int3{r_nc[0], r_cc[1], r_nc[2]} + grid.ldims + grid.ibn;

  lz = -Int3{l_nc[0], l_nc[1], l_cc[2]} + grid.ibn;
  rz = Int3{r_nc[0], r_nc[1], r_cc[2]} + grid.ldims + grid.ibn;
}

} // namespace detail

// ----------------------------------------------------------------------
// correct
//
// Do the modified marder correction (See eq.(5, 7, 9, 10) in Mardahl and
// Verboncoeur, CPC, 1997)

template <typename S>
inline void correct(const Grid_t& grid, S& mflds, const Int3& mflds_ib, S& mf,
                    const Int3& mf_ib, typename S::value_type diffusion)
{
  using real_t = typename S::value_type;
  using real3_t = gt::sarray<real_t, 3>;

  assert(mflds_ib == -grid.ibn);
  assert(mf_ib == -grid.ibn);

  real3_t fac = {.5f * grid.dt * diffusion / grid.domain.dx[0],
                 .5f * grid.dt * diffusion / grid.domain.dx[1],
                 .5f * grid.dt * diffusion / grid.domain.dx[2]};

  for (int p = 0; p < grid.n_patches(); p++) {
    Int3 lx, rx, ly, ry, lz, rz;
    detail::find_limits(grid, p, lx, rx, ly, ry, lz, rz);

    if (!grid.isInvar(0)) {
      Int3 l = lx, r = rx;
      auto ex = mflds.view(_all, _all, _all, EX, p);
      auto res = mf.view(_all, _all, _all, 0, p);
      ex.view(_s(l[0], r[0]), _s(l[1], r[1]), _s(l[2], r[2])) =
        ex.view(_s(l[0], r[0]), _s(l[1], r[1]), _s(l[2], r[2])) +
        (res.view(_s(l[0] + 1, r[0] + 1), _s(l[1], r[1]), _s(l[2], r[2])) -
         res.view(_s(l[0], r[0]), _s(l[1], r[1]), _s(l[2], r[2]))) *
          fac[0];
    }

    {
      Int3 l = ly, r = ry;
      auto ey = mflds.view(_all, _all, _all, EY, p);
      auto res = mf.view(_all, _all, _all, 0, p);
      ey.view(_s(l[0], r[0]), _s(l[1], r[1]), _s(l[2], r[2])) =
        ey.view(_s(l[0], r[0]), _s(l[1], r[1]), _s(l[2], r[2])) +
        (res.view(_s(l[0], r[0]), _s(l[1] + 1, r[1] + 1), _s(l[2], r[2])) -
         res.view(_s(l[0], r[0]), _s(l[1], r[1]), _s(l[2], r[2]))) *
          fac[1];
    }

    {
      Int3 l = lz, r = rz;
      auto ez = mflds.view(_all, _all, _all, EZ, p);
      auto res = mf.view(_all, _all, _all, 0, p);
      ez.view(_s(l[0], r[0]), _s(l[1], r[1]), _s(l[2], r[2])) =
        ez.view(_s(l[0], r[0]), _s(l[1], r[1]), _s(l[2], r[2])) +
        (res.view(_s(l[0], r[0]), _s(l[1], r[1]), _s(l[2] + 1, r[2] + 1)) -
         res.view(_s(l[0], r[0]), _s(l[1], r[1]), _s(l[2], r[2]))) *
          fac[2];
    }
  }
}

#ifdef USE_CUDA

inline void correct(const Grid_t& grid, MfieldsCuda::Storage& mflds,
                    const Int3& mflds_ib, MfieldsCuda::Storage& mf,
                    const Int3& mf_ib,
                    typename MfieldsCuda::Storage::value_type diffusion)
{
  Float3 fac;
  fac[0] = .5 * grid.dt * diffusion / grid.domain.dx[0];
  fac[1] = .5 * grid.dt * diffusion / grid.domain.dx[1];
  fac[2] = .5 * grid.dt * diffusion / grid.domain.dx[2];

  // OPT, do all patches in one kernel
  for (int p = 0; p < grid.n_patches(); p++) {
    assert(mflds.ib() == -grid.ibn());
    assert(mf.ib() == -grid.ibn());
    Int3 lx, rx, ly, ry, lz, rz;
    detail::find_limits(grid, p, lx, rx, ly, ry, lz, rz);

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
    psc::marder::correct(grid, mflds.storage(), mflds.ib(), res_.storage(),
                         res_.ib(), diffusion);
  }

  // ----------------------------------------------------------------------
  // operator()

  void operator()(MfieldsState& mflds, Mparticles& mprts)
  {
    static int pr;
    if (!pr) {
      pr = prof_register("marder", 1., 0, 0);
    }

    prof_start(pr);
    const auto& grid = mprts.grid();
    auto item_rho = Item_rho_t{grid};
    auto&& rho = psc::mflds::interior(grid, item_rho(mprts));

    // need to fill ghost cells first (should be unnecessary with only variant
    // 1) FIXME
    bnd_.fill_ghosts(mflds, EX, EX + 3);

    for (int i = 0; i < loop_; i++) {
      calc_aid_fields(mflds, rho);
      print_max(res_);
      correct(mflds);
      bnd_.fill_ghosts(mflds, EX, EX + 3);
    }
    prof_stop(pr);
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

template <typename MP, typename MFS, typename MF, typename D>
using Marder_ = MarderCommon<MP, MFS, MF, D,
                             Moment_rho_1st_nc<typename MF::Storage, D>, Bnd_>;
