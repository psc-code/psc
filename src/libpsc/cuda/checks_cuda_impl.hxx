
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

  ChecksCuda(const Grid_t& grid, MPI_Comm comm, const ChecksParams& params)
    : ChecksParams(params), continuity_{params}
  {}

  void continuity_before_particle_push(Mparticles& mprts)
  {
    continuity_.before_particle_push(mprts);
  }

  template <typename MfieldsState>
  void continuity_after_particle_push(Mparticles& mprts, MfieldsState& mflds)
  {
    continuity_.after_particle_push(mprts, mflds);
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
    auto rho = psc::mflds::interior(grid, item_rho(mprts));
    auto dive = psc::mflds::interior(grid, item_dive(mflds));
    auto&& h_rho = gt::host_mirror(rho);
    auto&& h_dive = gt::host_mirror(dive);
    gt::copy(gt::eval(rho), h_rho);
    gt::copy(gt::eval(dive), h_dive);

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
          double val_rho = h_rho(i, j, k, 0, p);
          double val_dive = h_dive(i, j, k, 0, p);
          max_err = std::max(max_err, std::abs(val_dive - val_rho));
          if (std::abs(val_dive - val_rho) > gauss_threshold) {
            printf("(%d,%d,%d): %g -- %g diff %g\n", i, j, k, val_dive, val_rho,
                   val_dive - val_rho);
          }
        }
      });
    }

    // find global max
    double tmp = max_err;
    MPI_Allreduce(&tmp, &max_err, 1, MPI_DOUBLE, MPI_MAX, grid.comm());

    if (gauss_verbose || max_err >= gauss_threshold) {
      mpi_printf(grid.comm(), "gauss: max_err = %g (thres %g)\n", max_err,
                 gauss_threshold);
    }

    if (gauss_dump_always || max_err >= gauss_threshold) {
      if (!writer_gauss_) {
        writer_gauss_.open("gauss");
      }
      writer_gauss_.begin_step(grid.timestep(), grid.timestep() * grid.dt);
      writer_gauss_.write(h_rho, grid, "rho", {"rho"});
      writer_gauss_.write(h_dive, grid, "dive", {"dive"});
      writer_gauss_.end_step();
    }

    assert(max_err < gauss_threshold);
  }

private:
  psc::checks::continuity<storage_type, Moment_t> continuity_;
  WriterDefault writer_gauss_;
};
