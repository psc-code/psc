
#pragma once

#include "checks.hxx"

#include "fields_item_moments_1st_cuda.hxx"

#ifdef PSC_HAVE_ADIOS2

#include "writer_adios2.hxx"
using WriterDefault = WriterADIOS2;

#else

#include "writer_mrc.hxx"
using WriterDefault = WriterMRC;

#endif

// FIXME!!! need to get Dim from Mparticles directly!

template <typename BS>
struct BS_to_Dim;

template <>
struct BS_to_Dim<BS144>
{
  using Dim = dim_yz;
};

template <>
struct BS_to_Dim<BS444>
{
  using Dim = dim_xyz;
};

template <typename Mparticles>
struct ChecksCuda
  : ChecksBase
  , ChecksParams
{
  using MfieldsState = MfieldsStateSingle;
  using Mfields = MfieldsSingle;
  using BS = typename Mparticles::BS;
  using Dim = typename BS_to_Dim<BS>::Dim;
  using Moment_t = Moment_rho_1st_nc_cuda<Dim>;
  using real_t = typename Mfields::real_t;

  ChecksCuda(const Grid_t& grid, MPI_Comm comm, const ChecksParams& params)
    : ChecksParams(params)
  {}

  void continuity_before_particle_push(Mparticles& mprts)
  {
    const auto& grid = mprts.grid();
    if (continuity_every_step <= 0 ||
        grid.timestep() % continuity_every_step != 0) {
      return;
    }

    Moment_t item_rho{grid};
    item_rho(mprts);
    rho_m_gt_ = psc::interior(item_rho.storage(), item_rho.ib());
  }

  void continuity_after_particle_push(Mparticles& mprts,
                                      MfieldsStateCuda& mflds)
  {
    const auto& grid = mprts.grid();
    if (continuity_every_step <= 0 ||
        grid.timestep() % continuity_every_step != 0) {
      return;
    }

    Moment_t item_rho{grid};
    item_rho(mprts);

    auto item_divj = Item_divj<MfieldsStateCuda>(mflds);

    auto d_rho_p = psc::interior(item_rho.storage(), item_rho.ib());
    auto&& rho_p = gt::host_mirror(d_rho_p);
    auto&& rho_m = gt::host_mirror(rho_m_gt_);
    auto&& h_divj = gt::host_mirror(item_divj.gt());
    gt::copy(gt::eval(d_rho_p), rho_p);
    gt::copy(gt::eval(rho_m_gt_), rho_m);
    gt::copy(gt::eval(item_divj.gt()), h_divj);

    auto&& d_rho = rho_p - rho_m;
    auto&& dt_divj = grid.dt * h_divj;

    double eps = continuity_threshold;
    double max_err = 0.;
    for (int p = 0; p < grid.n_patches(); p++) {
      grid.Foreach_3d(0, 0, [&](int i, int j, int k) {
        double _d_rho = d_rho(i, j, k, 0, p);
        double _dt_divj = dt_divj(i, j, k, 0, p);
        max_err = std::max(max_err, std::abs(_d_rho + _dt_divj));
        if (std::abs(_d_rho + _dt_divj) > eps) {
          mprintf("p%d (%d,%d,%d): %g -- %g diff %g\n", p, i, j, k, _d_rho,
                  -_dt_divj, _d_rho + _dt_divj);
        }
      });
    }

    // find global max
    double tmp = max_err;
    MPI_Allreduce(&tmp, &max_err, 1, MPI_DOUBLE, MPI_MAX, grid.comm());

    if (continuity_verbose || max_err >= eps) {
      mpi_printf(grid.comm(), "continuity: max_err = %g (thres %g)\n", max_err,
                 eps);
    }

    if (continuity_dump_always || max_err >= eps) {
      static WriterDefault writer;
      if (!writer) {
        writer.open("continuity");
      }
      writer.begin_step(grid.timestep(), grid.timestep() * grid.dt);
      writer.write(dt_divj, grid, "div_j", {"div_j"});
      writer.write(d_rho, grid, "d_rho", {"d_rho"});
      writer.end_step();
    }
    MPI_Barrier(grid.comm());

    assert(max_err < eps);
  }

  // ======================================================================
  // psc_checks: Gauss's Law

  // ----------------------------------------------------------------------
  // gauss

  void gauss(Mparticles& mprts, MfieldsStateCuda& mflds)
  {
    const auto& grid = mprts.grid();
    if (gauss_every_step <= 0 || grid.timestep() % gauss_every_step != 0) {
      return;
    }

    auto item_rho = Moment_t{grid};
    item_rho(mprts);
    auto item_dive = Item_dive<MfieldsStateCuda>(mflds);

    auto d_rho = psc::interior(item_rho.storage(), item_rho.ib());
    auto&& rho = gt::host_mirror(d_rho);
    auto&& dive = gt::host_mirror(item_dive.gt());
    gt::copy(gt::eval(d_rho), rho);
    gt::copy(gt::eval(item_dive.gt()), dive);

    double eps = gauss_threshold;
    double max_err = 0.;
    for (int p = 0; p < grid.n_patches(); p++) {
      int l[3] = {0, 0, 0}, r[3] = {0, 0, 0};
      for (int d = 0; d < 3; d++) {
        if (grid.bc.fld_lo[d] == BND_FLD_CONDUCTING_WALL &&
            grid.atBoundaryLo(p, d)) {
          l[d] = 1;
        }
      }

      grid.Foreach_3d(0, 0, [&](int jx, int jy, int jz) {
        if (jy < l[1] || jz < l[2] || jy >= grid.ldims[1] - r[1] ||
            jz >= grid.ldims[2] - r[2]) {
          // nothing
        } else {
          double v_rho = rho(jx, jy, jz, 0, p);
          double v_dive = dive(jx, jy, jz, 0, p);
          max_err = fmax(max_err, fabs(v_dive - v_rho));
#if 1
          if (fabs(v_dive - v_rho) > eps) {
            printf("(%d,%d,%d): %g -- %g diff %g\n", jx, jy, jz, v_dive, v_rho,
                   v_dive - v_rho);
          }
#endif
        }
      });
    }

    // find global max
    double tmp = max_err;
    MPI_Allreduce(&tmp, &max_err, 1, MPI_DOUBLE, MPI_MAX, grid.comm());

    if (gauss_verbose || max_err >= eps) {
      mpi_printf(grid.comm(), "gauss: max_err = %g (thres %g)\n", max_err, eps);
    }

    if (gauss_dump_always || max_err >= eps) {
      static WriterDefault writer;
      if (!writer) {
        writer.open("gauss");
      }
      writer.begin_step(grid.timestep(), grid.timestep() * grid.dt);
      writer.write(rho, grid, "rho", {"rho"});
      writer.write(dive, grid, "dive", {"dive"});
      writer.end_step();
    }

    MPI_Barrier(grid.comm());
    assert(max_err < eps);
  }

private:
  gt::gtensor<real_t, 5, gt::space::device> rho_m_gt_;
};
