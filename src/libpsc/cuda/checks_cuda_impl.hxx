
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

    double max_err = 0.;
    for (int p = 0; p < grid.n_patches(); p++) {
      int l[3] = {0, 0, 0}, r[3] = {0, 0, 0};
      for (int d = 0; d < 3; d++) {
        if (grid.bc.fld_lo[d] == BND_FLD_CONDUCTING_WALL &&
            grid.atBoundaryLo(p, d)) {
          l[d] = 1;
        }
      }

      auto patch_rho =
        rho.view(_s(l[0], -r[0]), _s(l[1], -r[1]), _s(l[2], -r[2]), 0, p);
      auto patch_dive =
        dive.view(_s(l[0], -r[0]), _s(l[1], -r[1]), _s(l[2], -r[2]), 0, p);

      auto patch_err = gt::norm_linf(patch_dive - patch_rho);
      max_err = std::max(max_err, patch_err);

      if (patch_err > gauss_threshold) {
        psc::helper::print_diff_3d(patch_rho, patch_dive, gauss_threshold);
      }
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
      writer_gauss_.write(rho, grid, "rho", {"rho"});
      writer_gauss_.write(dive, grid, "dive", {"dive"});
      writer_gauss_.end_step();
    }

    assert(max_err < gauss_threshold);
  }

private:
  psc::checks::continuity<storage_type, Moment_t> continuity_;
  WriterDefault writer_gauss_;
};
