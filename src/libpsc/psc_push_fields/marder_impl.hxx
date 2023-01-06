
#pragma once

#include "fields.hxx"
#include "writer_mrc.hxx"
#include "mpi_dtype_traits.hxx"
#include "../libpsc/psc_output_fields/fields_item_fields.hxx"
#include "../libpsc/psc_output_fields/fields_item_moments_1st.hxx"
#include "../libpsc/psc_bnd/psc_bnd_impl.hxx"

#include <gtensor/reductions.h>
#include <mrc_io.h>

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

template <typename E1, typename E2>
inline void correct(const Grid_t& grid, E1& efield, const Int3& efield_ib,
                    const E2& mf, const Int3& mf_ib,
                    typename E1::value_type diffusion)
{
  using real_t = gt::expr_value_type<E1>;
  using real3_t = gt::sarray<real_t, 3>;

  assert(efield_ib == -grid.ibn);
  assert(mf_ib == -grid.ibn);

  real3_t fac = {.5f * grid.dt * diffusion / grid.domain.dx[0],
                 .5f * grid.dt * diffusion / grid.domain.dx[1],
                 .5f * grid.dt * diffusion / grid.domain.dx[2]};

  for (int p = 0; p < grid.n_patches(); p++) {
    Int3 lx, rx, ly, ry, lz, rz;
    detail::find_limits(grid, p, lx, rx, ly, ry, lz, rz);

    if (!grid.isInvar(0)) {
      Int3 l = lx, r = rx;
      auto ex = efield.view(_all, _all, _all, 0, p);
      auto res = mf.view(_all, _all, _all, 0, p);
      ex.view(_s(l[0], r[0]), _s(l[1], r[1]), _s(l[2], r[2])) =
        ex.view(_s(l[0], r[0]), _s(l[1], r[1]), _s(l[2], r[2])) +
        (res.view(_s(l[0] + 1, r[0] + 1), _s(l[1], r[1]), _s(l[2], r[2])) -
         res.view(_s(l[0], r[0]), _s(l[1], r[1]), _s(l[2], r[2]))) *
          fac[0];
    }

    {
      Int3 l = ly, r = ry;
      auto ey = efield.view(_all, _all, _all, 1, p);
      auto res = mf.view(_all, _all, _all, 0, p);
      ey.view(_s(l[0], r[0]), _s(l[1], r[1]), _s(l[2], r[2])) =
        ey.view(_s(l[0], r[0]), _s(l[1], r[1]), _s(l[2], r[2])) +
        (res.view(_s(l[0], r[0]), _s(l[1] + 1, r[1] + 1), _s(l[2], r[2])) -
         res.view(_s(l[0], r[0]), _s(l[1], r[1]), _s(l[2], r[2]))) *
          fac[1];
    }

    {
      Int3 l = lz, r = rz;
      auto ez = efield.view(_all, _all, _all, 2, p);
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

template <typename E1, typename E2>
inline void cuda_marder_correct_yz(E1& efield, E2& res, Float3 fac, Int3 ly,
                                   Int3 ry, Int3 lz, Int3 rz)
{
  auto k_efield = efield.to_kernel();
  auto k_res = res.to_kernel();
  gt::launch<2>(
    {k_efield.shape(1), k_efield.shape(2)}, GT_LAMBDA(int iy, int iz) {
      if ((iy >= ly[1] && iy < ry[1]) && (iz >= ly[2] && iz < ry[2])) {
        k_efield(0, iy, iz, 1) =
          k_efield(0, iy, iz, 1) +
          fac[1] * (k_res(0, iy + 1, iz) - k_res(0, iy, iz));
      }

      if ((iy >= lz[1] && iy < rz[1]) && (iz >= lz[2] && iz < rz[2])) {
        k_efield(0, iy, iz, 2) =
          k_efield(0, iy, iz, 2) +
          fac[2] * (k_res(0, iy, iz + 1) - k_res(0, iy, iz));
      }
    });
  cuda_sync_if_enabled();
}

template <typename E1, typename E2>
inline void cuda_marder_correct_xyz(E1& efield, E2& res, Float3 fac, Int3 lx,
                                    Int3 rx, Int3 ly, Int3 ry, Int3 lz, Int3 rz)
{
  auto k_efield = efield.to_kernel();
  auto k_res = res.to_kernel();
  gt::launch<3>(
    {k_efield.shape(0), k_efield.shape(1), k_efield.shape(2)},
    GT_LAMBDA(int ix, int iy, int iz) {
      if ((ix >= lx[0] && ix < rx[0]) && (iy >= lx[1] && iy < rx[1]) &&
          (iz >= lx[2] && iz < rx[2])) {
        k_efield(ix, iy, iz, 0) =
          k_efield(ix, iy, iz, 0) +
          fac[0] * (k_res(ix, iy + 1, iz) - k_res(ix, iy, iz));
      }

      if ((ix >= ly[0] && ix < ry[0]) && (iy >= ly[1] && iy < ry[1]) &&
          (iz >= ly[2] && iz < ry[2])) {
        k_efield(ix, iy, iz, 1) =
          k_efield(ix, iy, iz, 1) +
          fac[1] * (k_res(ix, iy + 1, iz) - k_res(ix, iy, iz));
      }

      if ((ix >= lz[0] && ix < rz[0]) && (iy >= lz[1] && iy < rz[1]) &&
          (iz >= lz[2] && iz < rz[2])) {
        k_efield(ix, iy, iz, 2) =
          k_efield(ix, iy, iz, 2) +
          fac[2] * (k_res(ix, iy, iz + 1) - k_res(ix, iy, iz));
      }
    });
  cuda_sync_if_enabled();
}

template <typename E1>
inline void correct(const Grid_t& grid, E1& efield, const Int3& efield_ib,
                    MfieldsCuda::Storage& mf, const Int3& mf_ib,
                    typename MfieldsCuda::Storage::value_type diffusion)
{
  Float3 fac;
  fac[0] = .5 * grid.dt * diffusion / grid.domain.dx[0];
  fac[1] = .5 * grid.dt * diffusion / grid.domain.dx[1];
  fac[2] = .5 * grid.dt * diffusion / grid.domain.dx[2];

  assert(efield_ib == -grid.ibn);
  assert(mf_ib == -grid.ibn);
  // OPT, do all patches in one kernel
  for (int p = 0; p < grid.n_patches(); p++) {
    Int3 lx, rx, ly, ry, lz, rz;
    detail::find_limits(grid, p, lx, rx, ly, ry, lz, rz);

    auto p_efield = efield.view(_all, _all, _all, _all, p);
    auto p_res = mf.view(_all, _all, _all, 0, p);
    if (grid.isInvar(0)) {
      cuda_marder_correct_yz(p_efield, p_res, fac, ly, ry, lz, rz);
    } else {
      cuda_marder_correct_xyz(p_efield, p_res, fac, lx, rx, ly, ry, lz, rz);
    }
  }
}
#endif

} // namespace marder
} // namespace psc

template <typename S, typename D, typename ITEM_RHO, typename BND>
class MarderCommon
{
public:
  using storage_type = S;
  using dim_t = D;
  using Item_rho_t = ITEM_RHO;
  using Bnd = BND;
  using real_t = typename storage_type::value_type;

  // FIXME: checkpointing won't properly restore state

  MarderCommon(const Grid_t& grid, real_t diffusion, int loop, bool dump)
    : diffusion_{diffusion}, loop_{loop}, dump_{dump}, bnd_{grid, grid.ibn}
  {
    if (dump_) {
      io_.open("marder");
    }
  }

  // ----------------------------------------------------------------------
  // print_progress

  template <typename E1, typename E2, typename E3>
  void print_progress(const Grid_t& grid, const E1& rho, const E2& dive,
                      const E3& res)
  {
    real_t max_err = gt::norm_linf(res);
    MPI_Allreduce(MPI_IN_PLACE, &max_err, 1,
                  MpiDtypeTraits<gt::expr_value_type<E3>>::value(), MPI_MAX,
                  grid.comm());
    mpi_printf(grid.comm(), "marder: err %g\n", max_err);

    if (dump_) {
      static int cnt;
      io_.begin_step(cnt, cnt);
      cnt++;
      io_.write(rho, grid, "rho", {"rho"});
      io_.write(dive, grid, "dive", {"dive"});
      io_.end_step();
    }
  }

  // ----------------------------------------------------------------------
  // operator()
  //
  // Do the modified marder correction (See eq.(5, 7, 9, 10) in Mardahl and
  // Verboncoeur, CPC, 1997)

  template <typename Mparticles>
  void operator()(const Grid_t& grid, storage_type& mflds, const Int3& mflds_ib,
                  Mparticles& mprts)
  {
    auto efield = mflds.view(_all, _all, _all, _s(EX, EX + 3), _all);
    auto efield_ib = mflds_ib;

    double inv_sum = 0.;
    for (int d = 0; d < 3; d++) {
      if (!grid.isInvar(d)) {
        inv_sum += 1. / sqr(grid.domain.dx[d]);
      }
    }
    double diffusion_max = 1. / 2. / (.5 * grid.dt) / inv_sum;
    double diffusion = diffusion_max * diffusion_;

    auto item_rho = Item_rho_t{grid};
    auto rho = psc::mflds::interior(grid, item_rho(mprts));

    // FIXME: it's not well defined whether E ghost points are filled at entry,
    // and expected to be filled when done, so we're playing it safe for the
    // time being.
    for (int i = 0; i < loop_; i++) {
      bnd_.fill_ghosts(grid, mflds, mflds_ib, EX, EX + 3);
      auto dive = psc::mflds::interior(grid, psc::item::div_nc(grid, efield));

      Int3 res_ib = -grid.ibn;
      auto res = storage_type{psc::mflds::make_shape(grid, 1, res_ib)};
      psc::mflds::interior(grid, res) = dive - rho;
      bnd_.fill_ghosts(grid, res, res_ib, 0, 1);

      print_progress(grid, rho, dive, res);

      psc::marder::correct(grid, efield, efield_ib, res, res_ib, diffusion);
    }
    bnd_.fill_ghosts(grid, mflds, mflds_ib, EX, EX + 3);
  }

  template <typename MfieldsState, typename Mparticles>
  void operator()(MfieldsState& mflds, Mparticles& mprts)
  {
    static int pr;
    if (!pr) {
      pr = prof_register("marder", 1., 0, 0);
    }

    prof_start(pr);
    (*this)(mprts.grid(), mflds.storage(), mflds.ib(), mprts);
    prof_stop(pr);
  }

  // private:
  real_t diffusion_; //< diffusion coefficient for Marder correction
  int loop_;         //< execute this many relaxation steps in a loop
  bool dump_;        //< dump div_E, rho
  Bnd bnd_;
  WriterMRC io_; //< for debug dumping
};

template <typename S, typename D>
using Marder_ = MarderCommon<S, D, Moment_rho_1st_nc<S, D>, Bnd_>;

#ifdef USE_CUDA

#include "psc_particles_single.h"
#include "mparticles_cuda.hxx"
#include "fields_item_moments_1st_cuda.hxx"

template <typename D>
using MarderCuda = MarderCommon<MfieldsStateCuda::Storage, D,
                                Moment_rho_1st_nc_cuda<D>, BndCuda3>;
#endif
