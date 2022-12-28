
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

template <typename MP, typename D>
struct ChecksCuda : ChecksParams
{
  using Mparticles = MP;
  using dim_t = D;
  using storage_type = MfieldsCuda::Storage;
  using Moment_t = Moment_rho_1st_nc_cuda<D>;

  // ----------------------------------------------------------------------
  // ctor

  ChecksCuda(const Grid_t& grid, MPI_Comm comm, const ChecksParams& params)
    : ChecksParams(params)
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
    continuity_.before_particle_push(mprts);
  }

  // ----------------------------------------------------------------------
  // continuity_after_particle_push

  template <typename MfieldsState>
  void continuity_after_particle_push(Mparticles& mprts, MfieldsState& mflds)
  {
    const auto& grid = mprts.grid();
    if (continuity_every_step <= 0 ||
        grid.timestep() % continuity_every_step != 0) {
      return;
    }

    auto item_rho = Moment_t{grid};
    auto item_divj = Item_divj<MfieldsState>{};

    auto d_rho =
      psc::mflds::interior(grid, item_rho(mprts)) - continuity_.rho_m_;
    auto dt_divj = grid.dt * divj;

    double local_err = gt::norm_linf(d_rho + dt_divj);
    // find global max
    double max_err;
    MPI_Allreduce(&local_err, &max_err, 1, MPI_DOUBLE, MPI_MAX, grid.comm());

    if (max_err >= continuity_threshold) {
      auto&& h_d_rho = gt::host_mirror(d_rho);
      auto&& h_dt_divj = gt::host_mirror(dt_divj);
      gt::copy(gt::eval(d_rho), h_d_rho);
      gt::copy(gt::eval(dt_divj), h_dt_divj);
      for (int p = 0; p < grid.n_patches(); p++) {
        grid.Foreach_3d(0, 0, [&](int i, int j, int k) {
          double val_d_rho = h_d_rho(i, j, k, 0, p);
          double val_dt_divj = h_dt_divj(i, j, k, 0, p);
          if (std::abs(val_d_rho + val_dt_divj) > continuity_threshold) {
            mprintf("p%d (%d,%d,%d): %g -- %g diff %g\n", p, i, j, k, val_d_rho,
                    -val_dt_divj, val_d_rho + val_dt_divj);
          }
        });
      }
    }

    if (continuity_verbose || max_err >= continuity_threshold) {
      mpi_printf(grid.comm(), "continuity: max_err = %g (thres %g)\n", max_err,
                 continuity_threshold);
    }

    if (continuity_dump_always || max_err >= continuity_threshold) {
      if (!writer_continuity_) {
        writer_continuity_.open("continuity");
      }
      writer_continuity_.begin_step(grid.timestep(), grid.timestep() * grid.dt);
      writer_continuity_.write(d_rho, grid, "d_rho", {"d_rho"});
      writer_continuity_.write(dt_divj, grid, "div_j", {"div_j"});
      writer_continuity_.end_step();
    }
    MPI_Barrier(grid.comm());

    assert(max_err < eps);
  }

  // ======================================================================
  // psc_checks: Gauss's Law

  // ----------------------------------------------------------------------
  // gauss

  template <typename MfieldsState>
  void gauss(Mparticles& mprts, MfieldsState& mflds)
  {
    const auto& grid = mprts.grid();
    if (gauss_every_step <= 0 || grid.timestep() % gauss_every_step != 0) {
      return;
    }

    auto item_rho = Moment_t{grid};
    auto item_dive = Item_dive<MfieldsState>{};

    auto d_rho = psc::mflds::interior(grid, item_rho(mprts));
    auto d_dive = item_dive(mflds);
    auto&& rho = gt::host_mirror(d_rho);
    auto&& dive = gt::host_mirror(d_dive);
    gt::copy(d_rho, rho);
    gt::copy(d_dive, dive);

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
      if (!writer_gauss_) {
        writer_gauss_.open("gauss");
      }
      writer_gauss_.begin_step(grid.timestep(), grid.timestep() * grid.dt);
      writer_gauss_.write(rho, grid, "rho", {"rho"});
      writer_gauss_.write(dive, grid, "dive", {"dive"});
      writer_gauss_.end_step();
    }

    MPI_Barrier(grid.comm());
    assert(max_err < eps);
  }

private:
  psc::checks::continuity<storage_type, Moment_t> continuity_;
  WriterDefault writer_continuity_;
  WriterDefault writer_gauss_;
};
