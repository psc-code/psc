
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
  template <typename Mparticles, typename Mfields, typename D>
  using Moment_rho_nc = Moment_rho_1st_nc<Mparticles, Mfields, D>;
};

struct checks_order_2nd
{
  template <typename Mparticles, typename Mfields, typename D>
  using Moment_rho_nc = Moment_rho_2nd_nc<Mparticles, Mfields>;
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
  using Moment_t =
    typename ORDER::template Moment_rho_nc<Mparticles, Mfields, dim_t>;

  // ----------------------------------------------------------------------
  // ctor

  Checks_(const Grid_t& grid, MPI_Comm comm, const ChecksParams& params)
    : ChecksParams(params),
      comm_{comm},
      rho_{grid, 1, {}},
      rho_m_{grid, 1, {}},
      rho_p_{grid, 1, {}},
      divj_{grid, 1, {}}
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

    rho_m_.storage() = Moment_t{mprts}.gt();
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

    rho_p_.storage() = Moment_t{mprts}.gt();
    auto item_divj = Item_divj<MfieldsState>(mflds);

    auto& d_rho = rho_p_;
    d_rho.storage() = d_rho.storage() - rho_m_.storage();

    divj_.storage() = item_divj.gt();
    divj_.storage() = grid.dt * divj_.storage();

    double eps = continuity_threshold;
    double max_err = 0.;
    for (int p = 0; p < divj_.n_patches(); p++) {
      auto D_rho = make_Fields3d<dim_xyz>(d_rho[p]);
      auto Div_J = make_Fields3d<dim_xyz>(divj_[p]);
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
    MPI_Allreduce(&tmp, &max_err, 1, MPI_DOUBLE, MPI_MAX, comm_);

    if (continuity_verbose || max_err >= eps) {
      mpi_printf(comm_, "continuity: max_err = %g (thres %g)\n", max_err, eps);
    }

    if (continuity_dump_always || max_err >= eps) {
      if (!writer_continuity_) {
        writer_continuity_.open("continuity");
      }
      writer_continuity_.begin_step(grid.timestep(), grid.timestep() * grid.dt);
      writer_continuity_.write(divj_.gt(), grid, "div_j", {"div_j"});
      writer_continuity_.write(d_rho.gt(), grid, "d_rho", {"d_rho"});
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

    rho_.storage() = Moment_t{mprts}.gt();
    auto dive = Item_dive<MfieldsState>(mflds);
    auto&& dive_gt = view_interior(dive.gt(), dive.ibn());

    double eps = gauss_threshold;
    double max_err = 0.;
    for (int p = 0; p < dive.n_patches(); p++) {
      auto Rho = make_Fields3d<dim_xyz>(rho_[p]);

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
          double v_dive = dive_gt(jx, jy, jz, 0, p);
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
    MPI_Allreduce(&tmp, &max_err, 1, MPI_DOUBLE, MPI_MAX, comm_);

    if (gauss_verbose || max_err >= eps) {
      mpi_printf(comm_, "gauss: max_err = %g (thres %g)\n", max_err, eps);
    }

    if (gauss_dump_always || max_err >= eps) {
      if (!writer_gauss_) {
        writer_gauss_.open("gauss");
      }
      writer_gauss_.begin_step(grid.timestep(), grid.timestep() * grid.dt);
      writer_gauss_.write(rho_.gt(), grid, "rho", {"rho"});
      writer_gauss_.write(dive.gt(), dive.grid(), dive.name(),
                          dive.comp_names());
      writer_gauss_.end_step();
    }

    assert(max_err < eps);
  }

  // state
  MPI_Comm comm_;
  Mfields rho_p_;
  Mfields rho_m_;
  Mfields rho_;
  Mfields divj_;
  WriterDefault writer_continuity_;
  WriterDefault writer_gauss_;
};
