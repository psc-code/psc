
#pragma once

#include "fields.hxx"
#include "fields_item.hxx"
#include "checks.hxx"
#include "../libpsc/psc_output_fields/fields_item_fields.hxx"
#include "../libpsc/psc_output_fields/psc_output_fields_item_moments_1st_nc.cxx"
#include "../libpsc/psc_output_fields/psc_output_fields_item_moments_2nd_nc.cxx"

#include <mrc_io.h>

#ifdef PSC_HAVE_ADIOS2

#include "writer_adios2.hxx"
using WriterDefault = WriterADIOS2;

#else

#include "writer_mrc.hxx"
using WriterDefault = WriterMRC;

#endif

struct checks_order_1st
{
  template <typename Mfields, typename D>
  using Moment_rho_nc = Moment_rho_1st_nc<Mfields, D>;
};

struct checks_order_2nd
{
  template <typename Mfields, typename D>
  using Moment_rho_nc = Moment_rho_2nd_nc<Mfields, D>;
};

template <typename _Mparticles, typename _MfieldsState, typename _Mfields,
          typename ORDER, typename D>
struct Checks_
  : ChecksParams
  , ChecksBase
{
  using Mparticles = _Mparticles;
  using MfieldsState = _MfieldsState;
  using Mfields = _Mfields;
  using dim_t = D;
  using real_t = typename Mfields::real_t;
  using Moment_t = typename ORDER::template Moment_rho_nc<Mfields, dim_t>;

  // ----------------------------------------------------------------------
  // ctor

  Checks_(const Grid_t& grid, MPI_Comm comm, const ChecksParams& params)
    : ChecksParams(params), comm_{comm}
  {}

  // ======================================================================
  // psc_checks: Charge Continuity

  // ----------------------------------------------------------------------
  // continuity_before_particle_push

  void continuity_before_particle_push(Mparticles& mprts)
  {
    const auto& grid = mprts.grid();
    if (continuity_every_step <= 0 ||
        grid.timestep() % continuity_every_step != 0) {
      return;
    }

    auto rho = Moment_t{grid};
    rho_m_gt_ = psc::mflds::interior(grid, rho(mprts));
  }

  // ----------------------------------------------------------------------
  // continuity_after_particle_push

  void continuity_after_particle_push(Mparticles& mprts, MfieldsState& mflds)
  {
    const auto& grid = mprts.grid();
    if (continuity_every_step <= 0 ||
        grid.timestep() % continuity_every_step != 0) {
      return;
    }

    auto item_rho_p = Moment_t{grid};
    auto item_divj = Item_divj<MfieldsState>(mflds);

    auto d_rho_gt = psc::mflds::interior(grid, item_rho_p(mprts)) - rho_m_gt_;
    auto dt_divj_gt = grid.dt * item_divj.gt();

    double eps = continuity_threshold;
    double max_err = 0.;
    for (int p = 0; p < grid.n_patches(); p++) {
      grid.Foreach_3d(0, 0, [&](int i, int j, int k) {
        double val_d_rho = d_rho_gt(i, j, k, 0, p);
        double val_dt_divj = dt_divj_gt(i, j, k, 0, p);
        max_err = std::max(max_err, std::abs(val_d_rho + val_dt_divj));
        if (std::abs(val_d_rho + val_dt_divj) > eps) {
          mprintf("p%d (%d,%d,%d): %g -- %g diff %g\n", p, i, j, k, val_d_rho,
                  -val_dt_divj, val_d_rho + val_dt_divj);
        }
      });
    }

    // find global max
    double tmp = max_err;
    MPI_Allreduce(&tmp, &max_err, 1, MPI_DOUBLE, MPI_MAX, comm_);

    if (continuity_verbose || max_err >= eps) {
      mpi_printf(comm_, "continuity: max_err = %g (thres %g)\n", max_err, eps);
    }

    if (continuity_dump_always || max_err >= eps) {
      if (!writer_continuity_) {
        writer_continuity_.open("continuity");
      }
      writer_continuity_.begin_step(grid.timestep(), grid.timestep() * grid.dt);
      writer_continuity_.write(dt_divj_gt, grid, "dt_divj", {"dt_divj"});
      writer_continuity_.write(d_rho_gt, grid, "d_rho", {"d_rho"});
      writer_continuity_.end_step();
      MPI_Barrier(grid.comm());
    }

    assert(max_err < eps);
  }

  // ======================================================================
  // psc_checks: Gauss's Law

  // ----------------------------------------------------------------------
  // gauss

  void gauss(Mparticles& mprts, MfieldsState& mflds)
  {
    const auto& grid = mprts.grid();
    if (gauss_every_step <= 0 || grid.timestep() % gauss_every_step != 0) {
      return;
    }

    auto item_rho = Moment_t{grid};
    auto dive = Item_dive<MfieldsState>(mflds);
    auto rho_gt = psc::mflds::interior(grid, item_rho(mprts));
    auto dive_gt = psc::mflds::interior(grid, dive.storage());

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

      grid.Foreach_3d(0, 0, [&](int i, int j, int k) {
        if (j < l[1] || k < l[2] || j >= grid.ldims[1] - r[1] ||
            k >= grid.ldims[2] - r[2]) {
          // do nothing
        } else {
          double v_rho = rho_gt(i, j, k, 0, p);
          double v_dive = dive_gt(i, j, k, 0, p);
          max_err = std::max(max_err, std::abs(v_dive - v_rho));
          if (std::abs(v_dive - v_rho) > eps) {
            printf("(%d,%d,%d): %g -- %g diff %g\n", i, j, k, v_dive, v_rho,
                   v_dive - v_rho);
          }
        }
      });
    }

    // find global max
    double tmp = max_err;
    MPI_Allreduce(&tmp, &max_err, 1, MPI_DOUBLE, MPI_MAX, comm_);

    if (gauss_verbose || max_err >= eps) {
      mpi_printf(comm_, "gauss: max_err = %g (thres %g)\n", max_err, eps);
    }

    if (gauss_dump_always || max_err >= eps) {
      if (!writer_gauss_) {
        writer_gauss_.open("gauss");
      }
      writer_gauss_.begin_step(grid.timestep(), grid.timestep() * grid.dt);
      writer_gauss_.write(rho_gt, grid, "rho", {"rho"});
      writer_gauss_.write(dive_gt, grid, dive.name(), dive.comp_names());
      writer_gauss_.end_step();
    }

    assert(max_err < eps);
  }

  // state
  MPI_Comm comm_;
  gt::gtensor<real_t, 5> rho_m_gt_;
  WriterDefault writer_continuity_;
  WriterDefault writer_gauss_;
};
