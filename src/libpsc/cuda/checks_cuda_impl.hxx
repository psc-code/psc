
#pragma once

#include "checks.hxx"

#include "fields_item_moments_1st_cuda.hxx"

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
  using Moment_t = Moment_rho_1st_nc_cuda<Mparticles, Dim>;

  ChecksCuda(const Grid_t& grid, MPI_Comm comm, const ChecksParams& params)
    : ChecksParams(params),
      item_rho_{grid},
      item_rho_m_{grid},
      item_rho_p_{grid},
      divj_{grid, 1, grid.ibn}
  {}

  void continuity_before_particle_push(Mparticles& mprts)
  {
    const auto& grid = mprts.grid();
    if (continuity_every_step <= 0 ||
        grid.timestep() % continuity_every_step != 0) {
      return;
    }

    item_rho_m_(mprts);
  }

  void continuity_after_particle_push(Mparticles& mprts,
                                      MfieldsStateCuda& mflds)
  {
    const auto& grid = mprts.grid();
    if (continuity_every_step <= 0 ||
        grid.timestep() % continuity_every_step != 0) {
      return;
    }

    item_rho_p_(mprts);

    auto& h_mflds = mflds.get_as<MfieldsState>(0, mflds._n_comps());
    auto item_divj = Item_divj<MfieldsState>(h_mflds);

    auto& dev_rho_p = item_rho_p_.result();
    auto& dev_rho_m = item_rho_m_.result();
    auto& rho_p = dev_rho_p.template get_as<Mfields>(0, 1);
    auto& rho_m = dev_rho_m.template get_as<Mfields>(0, 1);

    auto& d_rho = rho_p;
    d_rho.axpy(-1., rho_m);

    divj_.assign(item_divj);
    divj_.scale(grid.dt);

    double eps = continuity_threshold;
    double max_err = 0.;
    for (int p = 0; p < divj_.n_patches(); p++) {
      auto D_rho = d_rho[p];
      auto Div_J = divj_[p];
      grid.Foreach_3d(0, 0, [&](int jx, int jy, int jz) {
        double d_rho = D_rho(0, jx, jy, jz);
        double div_j = Div_J(0, jx, jy, jz);
        max_err = fmax(max_err, fabs(d_rho + div_j));
        if (fabs(d_rho + div_j) > eps) {
          mprintf("p%d (%d,%d,%d): %g -- %g diff %g\n", p, jx, jy, jz, d_rho,
                  -div_j, d_rho + div_j);
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
      static struct mrc_io* io;
      if (!io) {
        io = mrc_io_create(grid.comm());
        mrc_io_set_name(io, "mrc_io_continuity");
        mrc_io_set_param_string(io, "basename", "continuity");
        mrc_io_set_from_options(io);
        mrc_io_setup(io);
        mrc_io_view(io);
      }
      mrc_io_open(io, "w", grid.timestep(), grid.timestep() * grid.dt);
      divj_.write_as_mrc_fld(io, "div_j", {"div_j"});
      d_rho.write_as_mrc_fld(io, "d_rho", {"d_rho"});
      mrc_io_close(io);
    }

    assert(max_err < eps);
    dev_rho_p.put_as(rho_p, 0, 0);
    dev_rho_m.put_as(rho_m, 0, 0);
    mflds.put_as(h_mflds, 0, 0);
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

    item_rho_(mprts);

    auto& h_mflds = mflds.get_as<MfieldsState>(0, mflds._n_comps());

    auto dive = Item_dive<MfieldsState>(h_mflds);
    auto& dev_rho = item_rho_.result();

    auto& rho = dev_rho.template get_as<Mfields>(0, 1);
    double eps = gauss_threshold;
    double max_err = 0.;
    for (int p = 0; p < dive.n_patches(); p++) {
      auto Rho = rho[p];

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
          double v_rho = Rho(0, jx, jy, jz);
          double v_dive = dive(0, {jx, jy, jz}, p);
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
      static struct mrc_io* io;
      if (!io) {
        io = mrc_io_create(grid.comm());
        mrc_io_set_name(io, "mrc_io_gauss");
        mrc_io_set_param_string(io, "basename", "gauss");
        mrc_io_set_from_options(io);
        mrc_io_setup(io);
        mrc_io_view(io);
      }
      mrc_io_open(io, "w", grid.timestep(), grid.timestep() * grid.dt);
      rho.write_as_mrc_fld(io, "rho", {"rho"});
      MrcIo::write_mflds(io, adaptMfields(dive), dive.grid(), dive.name(),
			 dive.comp_names());
      mrc_io_close(io);
    }

    assert(max_err < eps);
    dev_rho.put_as(rho, 0, 0);
    mflds.put_as(h_mflds, 0, 0);
  }

private:
  Moment_t item_rho_p_;
  Moment_t item_rho_m_;
  Moment_t item_rho_;
  Mfields divj_;
};
