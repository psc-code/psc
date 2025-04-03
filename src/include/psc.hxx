
#pragma once

#include <mrc_profile.h>
#include <DiagEnergies.h>

#include <particles.hxx>
#include <setup_particles.hxx>

#include "../libpsc/vpic/fields_item_vpic.hxx"
#include <checks_params.hxx>
#include <output_particles.hxx>
#include <push_particles.hxx>

#include "checkpoint.hxx"
#ifdef USE_CUDA
#include "../libpsc/cuda/mparticles_cuda.hxx"
#include "../libpsc/cuda/mparticles_cuda.inl"
#include "psc_fields_cuda.h"
#include "psc_fields_cuda.inl"
#endif

#include <fstream>

#ifdef VPIC
#include "../libpsc/vpic/vpic_iface.h"
#endif

#ifndef VPIC
struct MaterialList;
#endif

#ifdef VPIC

// FIXME, global variables are bad...

using VpicConfig = VpicConfigPsc;
using Mparticles = typename VpicConfig::Mparticles;
using MfieldsState = typename VpicConfig::MfieldsState;
using MfieldsHydro = MfieldsHydroQ<typename VpicConfig::Grid>;
using MfieldsInterpolator = typename VpicConfig::MfieldsInterpolator;
using MfieldsAccumulator = typename VpicConfig::MfieldsAccumulator;
using Grid = typename MfieldsState::Grid;
using ParticleBcList = typename Mparticles::ParticleBcList;
using MaterialList = typename MfieldsState::MaterialList;
using Material = typename MaterialList::Material;
using OutputHydro = OutputHydroQ<Mparticles, MfieldsHydro, MfieldsInterpolator>;
using DiagMixin =
  NoneDiagMixin<Mparticles, MfieldsState, MfieldsInterpolator, MfieldsHydro>;

Grid* vgrid;
std::unique_ptr<MfieldsHydro> hydro;
std::unique_ptr<MfieldsInterpolator> interpolator;
std::unique_ptr<MfieldsAccumulator> accumulator;
ParticleBcList particle_bc_list;
DiagMixin diag_mixin;

#endif

// ======================================================================

// FIXME, handle mem_fraction better
namespace detail
{
template <typename Mparticles, typename Enable = void>
struct mem_fraction
{
  static double get(const Mparticles& mprts) { return 0.; }
};

template <typename Mparticles>
struct mem_fraction<
  Mparticles,
  gt::meta::void_t<decltype(std::declval<Mparticles>().mem_fraction())>>
{
  static double get(const Mparticles& mprts)
  {
    double res = mprts.mem_fraction();
    MPI_Allreduce(MPI_IN_PLACE, &res, 1, MPI_DOUBLE, MPI_MAX,
                  mprts.grid().comm());
    return res;
  }
};

} // namespace detail

template <typename Mparticles>
double mem_fraction(const Mparticles& mprts)
{
  return detail::mem_fraction<Mparticles>::get(mprts);
}

// ======================================================================
// PscParams

struct PscParams
{
  double cfl = .75;            // CFL number used to determine time step
  int nmax;                    // Number of timesteps to run
  double wallclock_limit = 0.; // Maximum wallclock time to run
  int write_checkpoint_every_step = 0;

  bool detailed_profiling =
    false;              // output profiling info for each process separately
  int stats_every = 10; // output timing and other info every so many steps

  int balance_interval = 0;
  double balance_mem_fraction =
    1.; // balance if mprts.mem_fraction() exceeds this limit

  int sort_interval = 0;
  int marder_interval = 0;
};

// ----------------------------------------------------------------------
// courant_length

inline double courant_length(const Grid_t::Domain& domain)
{
  double inv_sum = 0.;
  for (int d = 0; d < 3; d++) {
    if (!domain.isInvar(d)) {
      inv_sum += 1. / sqr(domain.dx[d]);
    }
  }
  if (!inv_sum) { // simulation has 0 dimensions (happens in some test?)
    inv_sum = 1.;
  }
  return sqrt(1. / inv_sum);
}

// ======================================================================
// Psc

template <typename PscConfig, typename Diagnostics, typename InjectParticles,
          typename ExtCurrent>
struct Psc
{
  using Mparticles = typename PscConfig::Mparticles;
  using MfieldsState = typename PscConfig::MfieldsState;
  using Balance = typename PscConfig::Balance;
  using Sort = typename PscConfig::Sort;
  using Collision = typename PscConfig::Collision;
  using Checks = typename PscConfig::Checks;
  using Marder = typename PscConfig::Marder;
  using PushParticles = typename PscConfig::PushParticles;
  using PushFields = typename PscConfig::PushFields;
  using Bnd = typename PscConfig::Bnd;
  using BndFields = typename PscConfig::BndFields;
  using BndParticles = typename PscConfig::BndParticles;
  using Dim = typename PscConfig::Dim;

#ifdef VPIC
  using AccumulateOps = typename PushParticles::AccumulateOps;
#endif

  // ----------------------------------------------------------------------
  // ctor

  Psc(const PscParams& params, Grid_t& grid, MfieldsState& mflds,
      Mparticles& mprts, Balance& balance, Collision& collision, Checks& checks,
      Marder& marder, Diagnostics& diagnostics,
      InjectParticles& inject_particles, ExtCurrent& ext_current)
    : p_{params},
      grid_{&grid},
      mflds_{mflds},
      mprts_{mprts},
      balance_{balance},
      collision_{collision},
      checks_{checks},
      marder_{marder},
      bndp_{grid},
      diagnostics_{diagnostics},
      inject_particles_{inject_particles},
      ext_current_{ext_current},
      checkpointing_{params.write_checkpoint_every_step}
  {
    time_start_ = MPI_Wtime();

    assert(grid.isInvar(0) == Dim::InvarX::value);
    assert(grid.isInvar(1) == Dim::InvarY::value);
    assert(grid.isInvar(2) == Dim::InvarZ::value);

    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    log_.open("mem-" + std::to_string(rank) + ".log");

#ifdef USE_CUDA
    mem_stats_csv_header(log_);
#endif

    initialize_stats();
    initialize();
  }

  // ----------------------------------------------------------------------
  // initialize_stats

  void initialize_stats()
  {
    st_nr_particles = psc_stats_register("nr particles");
    st_time_step = psc_stats_register("time entire step");

    // generic stats categories
    st_time_particle = psc_stats_register("time particle update");
    st_time_field = psc_stats_register("time field update");
    st_time_comm = psc_stats_register("time communication");
    st_time_output = psc_stats_register("time output");

    // FIXME not quite the right place
    pr_time_step_no_comm = prof_register("time step w/o comm", 1., 0, 0);
  }

  // ----------------------------------------------------------------------
  // initialize

  void initialize()
  {
#ifdef VPIC
    initialize_vpic();
#else
    initialize_default();
#endif
    bndf_.fill_ghosts_H(mflds_);
    bnd_.fill_ghosts(mflds_, HX, HX + 3);

    bnd_.fill_ghosts(mflds_, JXI, JXI + 3);

    bndf_.fill_ghosts_E(mflds_);
    bnd_.fill_ghosts(mflds_, EX, EX + 3);

    // initial output / stats
    mpi_printf(grid().comm(), "Performing initial diagnostics.\n");
    diagnostics();

    psc_stats_val[st_nr_particles] = mprts_.size();

    print_status();

    mpi_printf(grid().comm(), "Initialization complete.\n");
  }

  // ----------------------------------------------------------------------
  // integrate

  void integrate()
  {
    static int pr;
    if (!pr) {
      pr = prof_register("psc_step", 1., 0, 0);
    }

    mpi_printf(grid().comm(), "*** Advancing\n");
    double elapsed = MPI_Wtime();

    while (grid().timestep() < p_.nmax) {
      prof_start(pr);
      psc_stats_start(st_time_step);

      checkpointing_(grid(), mprts_, mflds_);

      mpi_printf(grid().comm(),
                 "**** Step %d / %d, Code Time %g, Wall Time %g\n",
                 grid().timestep() + 1, p_.nmax, grid().timestep() * grid().dt,
                 MPI_Wtime() - time_start_);

      // prof_start(pr_time_step_no_comm);
      // prof_stop(
      //   pr_time_step_no_comm); // actual measurements are done w/ restart

      step();
      grid_->timestep_++; // FIXME, too hacky
#ifdef VPIC
      vgrid->step++;
      assert(vgrid->step == grid().timestep());
#endif

      diagnostics();

      psc_stats_stop(st_time_step);
      prof_stop(pr);

      psc_stats_val[st_nr_particles] = mprts_.size();

      if (grid().timestep() % p_.stats_every == 0) {
        print_status();
      }

      if (p_.wallclock_limit > 0.) {
        double wallclock_elapsed = MPI_Wtime() - time_start_;
        double wallclock_elapsed_max;
        MPI_Allreduce(&wallclock_elapsed, &wallclock_elapsed_max, 1, MPI_DOUBLE,
                      MPI_MAX, MPI_COMM_WORLD);

        if (wallclock_elapsed_max > p_.wallclock_limit) {
          mpi_printf(MPI_COMM_WORLD, "WARNING: Max wallclock time elapsed!\n");
          break;
        }
      }
    }

    checkpointing_.final(grid(), mprts_, mflds_);

    // FIXME, merge with existing handling of wallclock time
    elapsed = MPI_Wtime() - elapsed;

    int s = (int)elapsed, m = s / 60, h = m / 60, d = h / 24, w = d / 7;
    /**/ s -= m * 60, m -= h * 60, h -= d * 24, d -= w * 7;
    mpi_printf(grid().comm(),
               "*** Finished (%gs / %iw:%id:%ih:%im:%is elapsed)\n", elapsed, w,
               d, h, m, s);
  }

#ifdef VPIC
  // ----------------------------------------------------------------------
  // step_vpic

  void step_vpic()
  {
    static int pr_sort, pr_collision, pr_checks, pr_push_prts, pr_push_flds,
      pr_bndp, pr_bndf, pr_marder, pr_inject, pr_heating;
    if (!pr_sort) {
      pr_sort = prof_register("step_sort", 1., 0, 0);
      pr_collision = prof_register("step_collision", 1., 0, 0);
      pr_push_prts = prof_register("step_push_prts", 1., 0, 0);
      pr_push_flds = prof_register("step_push_flds", 1., 0, 0);
      pr_bndp = prof_register("step_bnd_prts", 1., 0, 0);
      pr_bndf = prof_register("step_bnd_flds", 1., 0, 0);
      pr_checks = prof_register("step_checks", 1., 0, 0);
      pr_marder = prof_register("step_marder", 1., 0, 0);
      pr_inject = prof_register("step_inject", 1., 0, 0);
      pr_heating = prof_register("step_heating", 1., 0, 0);
    }

    MPI_Comm comm = grid().comm();

    // x^{n+1/2}, p^{n}, E^{n+1/2}, B^{n+1/2}

    int timestep = grid().timestep();

    if (p_.balance_interval > 0 && timestep % p_.balance_interval == 0) {
      balance_(grid_, mprts_);
    }

    // prof_start(pr_time_step_no_comm);
    // prof_stop(pr_time_step_no_comm); // actual measurements are done w/
    // restart

    if (p_.sort_interval > 0 && timestep % p_.sort_interval == 0) {
      // mpi_printf(comm, "***** Sorting...\n");
      prof_start(pr_sort);
      sort_(mprts_);
      prof_stop(pr_sort);
    }

    if (collision_.interval() > 0 && timestep % collision_.interval() == 0) {
      mpi_printf(comm, "***** Performing collisions...\n");
      prof_start(pr_collision);
      collision_(mprts_);
      prof_stop(pr_collision);
    }

    // psc_bnd_particles_open_calc_moments(psc_->bnd_particles,
    // psc_->particles);

    checks_.continuity_before_particle_push(mprts_);

    // === particle propagation p^{n} -> p^{n+1}, x^{n+1/2} -> x^{n+3/2}
    prof_start(pr_push_prts);
    pushp_.push_mprts(mprts_, mflds_, *interpolator, *accumulator,
                      particle_bc_list, num_comm_round);
    prof_stop(pr_push_prts);
    // state is now: x^{n+3/2}, p^{n+1}, E^{n+1/2}, B^{n+1/2}, j^{n+1}

    // field propagation B^{n+1/2} -> B^{n+1}
    pushf_.push_H(mflds_, .5);
    // x^{n+3/2}, p^{n+1}, E^{n+1/2}, B^{n+1}, j^{n+1}

    prof_start(pr_bndp);
    bndp_(mprts_);
    prof_stop(pr_bndp);

    // field propagation E^{n+1/2} -> E^{n+3/2}

    // fill ghosts for H
    bndf_.fill_ghosts_H(mflds_);
    bnd_.fill_ghosts(mflds_, HX, HX + 3);

    // add and fill ghost for J
    bndf_.add_ghosts_J(mflds_);
    bnd_.add_ghosts(mflds_, JXI, JXI + 3);
    bnd_.fill_ghosts(mflds_, JXI, JXI + 3);

    // push E
    pushf_.push_E(mflds_, 1.);

    bndf_.fill_ghosts_E(mflds_);
    // if (pushf_->variant == 0) {
    bnd_.fill_ghosts(mflds_, EX, EX + 3);
    //}
    // x^{n+3/2}, p^{n+1}, E^{n+3/2}, B^{n+1}

    // field propagation B^{n+1} -> B^{n+3/2}
    // if (pushf_->variant == 0) {
    bndf_.fill_ghosts_E(mflds_);
    bnd_.fill_ghosts(mflds_, EX, EX + 3);
    //    }

    // push H
    pushf_.push_H(mflds_, .5);

    bndf_.fill_ghosts_H(mflds_);
    // if (pushf_->variant == 0) {
    bnd_.fill_ghosts(mflds_, HX, HX + 3);
    //}
    // x^{n+3/2}, p^{n+1}, E^{n+3/2}, B^{n+3/2}

    checks_.continuity_after_particle_push(mprts_, mflds_);

    // E at t^{n+3/2}, particles at t^{n+3/2}
    // B at t^{n+3/2} (Note: that is not it's natural time,
    // but div B should be == 0 at any time...)
    if (p_.marder_interval > 0 && timestep % p_.marder_interval == 0) {
      // mpi_printf(comm, "***** Performing Marder correction...\n");
      prof_start(pr_marder);
      marder_(mflds_, mprts_);
      prof_stop(pr_marder);
    }

    checks_.gauss(mprts_, mflds_);

#ifdef VPIC
    pushp_.load_interpolator(mprts_, mflds_, *interpolator);
#endif
  }
#endif

  // ----------------------------------------------------------------------
  // step_psc

  void step_psc()
  {
    using Dim = typename PscConfig::Dim;

    static int pr_sort, pr_collision, pr_checks, pr_push_prts, pr_push_flds,
      pr_bndp, pr_bndf, pr_marder, pr_inject_prts;
    if (!pr_sort) {
      pr_sort = prof_register("step_sort", 1., 0, 0);
      pr_collision = prof_register("step_collision", 1., 0, 0);
      pr_push_prts = prof_register("step_push_prts", 1., 0, 0);
      pr_inject_prts = prof_register("step_inject_prts", 1., 0, 0);
      pr_push_flds = prof_register("step_push_flds", 1., 0, 0);
      pr_bndp = prof_register("step_bnd_prts", 1., 0, 0);
      pr_bndf = prof_register("step_bnd_flds", 1., 0, 0);
      pr_checks = prof_register("step_checks", 1., 0, 0);
      pr_marder = prof_register("step_marder", 1., 0, 0);
    }

    // state is at: x^{n+1/2}, p^{n}, E^{n+1/2}, B^{n+1/2}
    MPI_Comm comm = grid().comm();
    int timestep = grid().timestep();

#ifdef USE_CUDA
    mem_stats_csv(log_, timestep, grid().n_patches(), mprts_.size());
#endif

    double mem_fraction = ::mem_fraction(mprts_);
    if (p_.balance_interval > 0 && (timestep % p_.balance_interval == 0 ||
                                    mem_fraction > p_.balance_mem_fraction)) {
      balance_(grid_, mprts_);
    }

    if (p_.sort_interval > 0 && timestep % p_.sort_interval == 0) {
      mpi_printf(comm, "***** Sorting...\n");
      prof_start(pr_sort);
      sort_(mprts_);
      prof_stop(pr_sort);
    }

    if (collision_.interval() > 0 && timestep % collision_.interval() == 0) {
      mpi_printf(comm, "***** Performing collisions...\n");
      prof_start(pr_collision);
      collision_(mprts_);
      prof_stop(pr_collision);
    }

    // === particle injection
    prof_start(pr_inject_prts);
    inject_particles();
    prof_stop(pr_inject_prts);

    if (checks_.continuity.every_step > 0 &&
        timestep % checks_.continuity.every_step == 0) {
      mpi_printf(comm, "***** Checking continuity...\n");
      prof_start(pr_checks);
      checks_.continuity_before_particle_push(mprts_);
      prof_stop(pr_checks);
    }

    // === particle propagation p^{n} -> p^{n+1}, x^{n+1/2} -> x^{n+3/2}
    mpi_printf(comm, "***** Pushing particles...\n");
    prof_start(pr_push_prts);
    pushp_.push_mprts(mprts_, mflds_);
    prof_stop(pr_push_prts);
    // state is now: x^{n+3/2}, p^{n+1}, E^{n+1/2}, B^{n+1/2}, j^{n+1}

    // === field propagation B^{n+1/2} -> B^{n+1}
    mpi_printf(comm, "***** Pushing B...\n");
    prof_start(pr_push_flds);
    pushf_.push_H(mflds_, .5, Dim{});
    prof_stop(pr_push_flds);
    // state is now: x^{n+3/2}, p^{n+1}, E^{n+1/2}, B^{n+1}, j^{n+1}

    mpi_printf(comm, "***** Bnd particles...\n");
    prof_start(pr_bndp);
    bndp_(mprts_);
    prof_stop(pr_bndp);

    // === field propagation E^{n+1/2} -> E^{n+3/2}
    mpi_printf(comm, "***** Push fields E\n");
    prof_start(pr_bndf);
#if 1
    bndf_.fill_ghosts_H(mflds_);
    bnd_.fill_ghosts(mflds_, HX, HX + 3);
#endif

    // === external current
    this->ext_current_(grid(), mflds_);

    bndf_.add_ghosts_J(mflds_);
    bnd_.add_ghosts(mflds_, JXI, JXI + 3);
    bnd_.fill_ghosts(mflds_, JXI, JXI + 3);
    prof_stop(pr_bndf);

    prof_restart(pr_push_flds);
    pushf_.push_E(mflds_, 1., Dim{});
    prof_stop(pr_push_flds);

#if 1
    prof_restart(pr_bndf);
    bndf_.fill_ghosts_E(mflds_);
    bnd_.fill_ghosts(mflds_, EX, EX + 3);
    prof_stop(pr_bndf);
#endif
    // state is now: x^{n+3/2}, p^{n+1}, E^{n+3/2}, B^{n+1}

    // === field propagation B^{n+1} -> B^{n+3/2}
    mpi_printf(comm, "***** Push fields B\n");
    prof_restart(pr_push_flds);
    pushf_.push_H(mflds_, .5, Dim{});
    prof_stop(pr_push_flds);

#if 1
    prof_start(pr_bndf);
    bndf_.fill_ghosts_H(mflds_);
    bnd_.fill_ghosts(mflds_, HX, HX + 3);
    prof_stop(pr_bndf);
    // state is now: x^{n+3/2}, p^{n+1}, E^{n+3/2}, B^{n+3/2}
#endif

    if (checks_.continuity.every_step > 0 &&
        timestep % checks_.continuity.every_step == 0) {
      prof_restart(pr_checks);
      checks_.continuity_after_particle_push(mprts_, mflds_);
      prof_stop(pr_checks);
    }

    // E at t^{n+3/2}, particles at t^{n+3/2}
    // B at t^{n+3/2} (Note: that is not its natural time,
    // but div B should be == 0 at any time...)
    if (p_.marder_interval > 0 && timestep % p_.marder_interval == 0) {
      mpi_printf(comm, "***** Performing Marder correction...\n");
      prof_start(pr_marder);
      marder_(mflds_, mprts_);
      prof_stop(pr_marder);
    }

    if (checks_.gauss.every_step > 0 &&
        timestep % checks_.gauss.every_step == 0) {
      prof_restart(pr_checks);
      checks_.gauss(mprts_, mflds_);
      prof_stop(pr_checks);
    }

    // psc_push_particles_prep(psc->push_particles, psc->particles, psc->flds);
  }

  void step()
  {
#ifdef VPIC
    step_vpic();
#else
    step_psc();
#endif
  }

  // ----------------------------------------------------------------------
  // inject_particles

  void inject_particles() { return this->inject_particles_(grid(), mprts_); }

private:
  // ----------------------------------------------------------------------
  // print_profiling

  void print_profiling()
  {
    int size;
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    if (!p_.detailed_profiling) {
      prof_print_mpi(MPI_COMM_WORLD);
    } else {
      int rank;
      MPI_Comm_rank(MPI_COMM_WORLD, &rank);
      for (int i = 0; i < size; i++) {
        if (i == rank) {
          mprintf("profile\n");
          prof_print();
        }
        MPI_Barrier(MPI_COMM_WORLD);
      }
    }
  }

#ifndef VPIC
  // ----------------------------------------------------------------------
  // initialize_default

  void initialize_default()
  {
    // pushp_.stagger(mprts, mflds); FIXME, vpic does it

    // checks_.gauss(mprts_, mflds_);
  }
#endif

#ifdef VPIC

  // ----------------------------------------------------------------------
  // initialize_vpic

  void initialize_vpic()
  {
    MPI_Comm comm = grid().comm();

    // Do some consistency checks on user initialized fields

    mpi_printf(comm, "Checking interdomain synchronization\n");
    double err = marder_.synchronize_tang_e_norm_b(mflds_);
    mpi_printf(comm, "Error = %g (arb units)\n", err);

    mpi_printf(comm, "Checking magnetic field divergence\n");
    marder_.compute_div_b_err(mflds_);
    err = marder_.compute_rms_div_b_err(mflds_);
    mpi_printf(comm, "RMS error = %e (charge/volume)\n", err);
    marder_.clean_div_b(mflds_);

    // Load fields not initialized by the user

    mpi_printf(comm, "Initializing radiation damping fields\n");
    TIC AccumulateOps::compute_curl_b(mflds_);
    TOC(compute_curl_b, 1);

    mpi_printf(comm, "Initializing bound charge density\n");
    marder_.clear_rhof(mflds_);
    marder_.accumulate_rho_p(mprts_, mflds_);
    marder_.synchronize_rho(mflds_);
    TIC AccumulateOps::compute_rhob(mflds_);
    TOC(compute_rhob, 1);

    // Internal sanity checks

    mpi_printf(comm, "Checking electric field divergence\n");
    marder_.compute_div_e_err(mflds_);
    err = marder_.compute_rms_div_e_err(mflds_);
    mpi_printf(comm, "RMS error = %e (charge/volume)\n", err);
    marder_.clean_div_e(mflds_);

    mpi_printf(comm, "Rechecking interdomain synchronization\n");
    err = marder_.synchronize_tang_e_norm_b(mflds_);
    mpi_printf(comm, "Error = %e (arb units)\n", err);

    mpi_printf(comm, "Uncentering particles\n");
    if (!mprts_.empty()) {
      pushp_.load_interpolator(mprts_, mflds_, *interpolator);
      pushp_.uncenter(mprts_, *interpolator);
    }
  }

#endif

  // ----------------------------------------------------------------------
  // diagnostics

  void diagnostics() { diagnostics_(mprts_, mflds_); }

  // ----------------------------------------------------------------------
  // print_status

  void print_status()
  {
#ifdef VPIC
#ifdef USE_VPIC
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    update_profile(rank == 0);
#endif
#endif
    psc_stats_log(grid().timestep());
    print_profiling();
  }

public:
  const Grid_t& grid() { return *grid_; }

private:
  double time_start_;
  PscParams p_;

protected:
  Grid_t* grid_;

  MfieldsState& mflds_;
  Mparticles& mprts_;

  Balance& balance_;
  Collision& collision_;
  Checks& checks_;
  Marder& marder_;
  Diagnostics& diagnostics_;
  InjectParticles& inject_particles_;
  ExtCurrent& ext_current_;

  Sort sort_;
  PushParticles pushp_;
  PushFields pushf_;
  Bnd bnd_;
  BndFields bndf_;
  BndParticles bndp_;

  Checkpointing checkpointing_;
  std::ofstream log_;

  // FIXME, maybe should be private
  // need to make sure derived class sets these (? -- or just leave them off by
  // default)
  int num_comm_round = {3};

  int st_nr_particles;
  int st_time_step;
};

// ======================================================================
// InjectParticlesNone

class InjectParticlesNone
{
public:
  template <typename Mparticles>
  void operator()(const Grid_t& grid, Mparticles& mprts)
  {}
};

namespace
{

InjectParticlesNone injectParticlesNone;

} // namespace

// ======================================================================
// ExtCurrentNone

class ExtCurrentNone
{
public:
  template <typename MfieldsState>
  void operator()(const Grid_t& grid, MfieldsState& mflds)
  {}
};

namespace
{

ExtCurrentNone extCurrentNone;

} // namespace
// ======================================================================
// makePscIntegrator

template <typename PscConfig, typename MfieldsState, typename Mparticles,
          typename Balance, typename Collision, typename Checks,
          typename Marder, typename Diagnostics,
          typename InjectParticles = InjectParticlesNone,
          typename ExtCurrent = ExtCurrentNone>
Psc<PscConfig, Diagnostics, InjectParticles, ExtCurrent> makePscIntegrator(
  const PscParams& params, Grid_t& grid, MfieldsState& mflds, Mparticles& mprts,
  Balance& balance, Collision& collision, Checks& checks, Marder& marder,
  Diagnostics& diagnostics,
  InjectParticles& inject_particles = injectParticlesNone,
  ExtCurrent& ext_current = extCurrentNone)
{
  return {params,     grid,   mflds,  mprts,       balance,
          collision,  checks, marder, diagnostics, inject_particles,
          ext_current};
}

// ======================================================================
// partitionAndSetupParticles

template <typename SetupParticles, typename Balance, typename Mparticles,
          typename InitFunc>
void partitionAndSetupParticles(SetupParticles& setup_particles,
                                Balance& balance, Grid_t*& grid_ptr,
                                Mparticles& mprts, InitFunc init_np)
{
  partitionParticles(setup_particles, balance, grid_ptr, mprts, init_np);
  setupParticles(setup_particles, mprts, init_np);
}

// ----------------------------------------------------------------------
// partitionParticles

template <typename SetupParticles, typename Balance, typename Mparticles,
          typename InitFunc>
void partitionParticles(SetupParticles& setup_particles, Balance& balance,
                        Grid_t*& grid_ptr, Mparticles& mprts, InitFunc init_np)
{
  auto comm = grid_ptr->comm();
  mpi_printf(comm, "**** Partitioning...\n");

  auto n_prts_by_patch = setup_particles.partition(*grid_ptr, init_np);

  balance.initial(grid_ptr, n_prts_by_patch);
  // !!! FIXME! grid is now invalid
  // balance::initial does not rebalance particles, because the old way of
  // doing this does't even have the particle data structure created yet --
  // FIXME?
  mprts.reset(*grid_ptr);
}

// ----------------------------------------------------------------------
// setupParticles

template <typename SetupParticles, typename Mparticles, typename InitFunc>
void setupParticles(SetupParticles& setup_particles, Mparticles& mprts,
                    InitFunc init_np)
{
  mpi_printf(MPI_COMM_WORLD, "**** Setting up particles...\n");
  setup_particles.setupParticles(mprts, init_np);
}

// ======================================================================
// VPIC-like stuff

// ----------------------------------------------------------------------
// define_periodic_grid

void define_periodic_grid(const double xl[3], const double xh[3],
                          const int gdims[3], const int np[3])
{
#ifdef VPIC
  // SimulationMixin::setTopology(np[0], np[1], np[2]); FIXME, needed for
  // vpic_simulation, I believe only because this info is written out in
  // diagnostics_run
  vgrid->partition_periodic_box(xl, xh, gdims, Int3::fromPointer(np));
#endif
}

// ----------------------------------------------------------------------
// set_domain_field_bc

void set_domain_field_bc(Int3 bnd, int bc)
{
#ifdef VPIC
  int boundary = BOUNDARY(bnd[0], bnd[1], bnd[2]);
  int fbc;
  switch (bc) {
    case BND_FLD_CONDUCTING_WALL: fbc = Grid::pec_fields; break;
    case BND_FLD_ABSORBING: fbc = Grid::absorb_fields; break;
    default: assert(0);
  }
  vgrid->set_fbc(boundary, fbc);
#endif
}

// ----------------------------------------------------------------------
// set_domain_particle_bc

void set_domain_particle_bc(Int3 bnd, int bc)
{
#ifdef VPIC
  int boundary = BOUNDARY(bnd[0], bnd[1], bnd[2]);
  int pbc;
  switch (bc) {
    case BND_PRT_REFLECTING: pbc = Grid::reflect_particles; break;
    case BND_PRT_ABSORBING: pbc = Grid::absorb_particles; break;
    default: assert(0);
  }
  vgrid->set_pbc(boundary, pbc);
#endif
}

void grid_setup_communication()
{
#ifdef VPIC
  assert(vgrid->nx && vgrid->ny && vgrid->ny);

  // Pre-size communications buffers. This is done to get most memory
  // allocation over with before the simulation starts running
  // FIXME, this isn't a great place. First, we shouldn't call mp
  // functions (semi-)directly. 2nd, whether we need these buffers depends
  // on b.c., which aren't yet known.

  // FIXME, this really isn't a good place to do this, as it requires layer
  // breaking knowledge of which communication will need the largest
  // buffers...
  int nx1 = vgrid->nx + 1, ny1 = vgrid->ny + 1, nz1 = vgrid->nz + 1;
  vgrid->mp_size_recv_buffer(
    BOUNDARY(-1, 0, 0), ny1 * nz1 * sizeof(typename MfieldsHydro::Element));
  vgrid->mp_size_recv_buffer(
    BOUNDARY(1, 0, 0), ny1 * nz1 * sizeof(typename MfieldsHydro::Element));
  vgrid->mp_size_recv_buffer(
    BOUNDARY(0, -1, 0), nz1 * nx1 * sizeof(typename MfieldsHydro::Element));
  vgrid->mp_size_recv_buffer(
    BOUNDARY(0, 1, 0), nz1 * nx1 * sizeof(typename MfieldsHydro::Element));
  vgrid->mp_size_recv_buffer(
    BOUNDARY(0, 0, -1), nx1 * ny1 * sizeof(typename MfieldsHydro::Element));
  vgrid->mp_size_recv_buffer(
    BOUNDARY(0, 0, 1), nx1 * ny1 * sizeof(typename MfieldsHydro::Element));

  vgrid->mp_size_send_buffer(
    BOUNDARY(-1, 0, 0), ny1 * nz1 * sizeof(typename MfieldsHydro::Element));
  vgrid->mp_size_send_buffer(
    BOUNDARY(1, 0, 0), ny1 * nz1 * sizeof(typename MfieldsHydro::Element));
  vgrid->mp_size_send_buffer(
    BOUNDARY(0, -1, 0), nz1 * nx1 * sizeof(typename MfieldsHydro::Element));
  vgrid->mp_size_send_buffer(
    BOUNDARY(0, 1, 0), nz1 * nx1 * sizeof(typename MfieldsHydro::Element));
  vgrid->mp_size_send_buffer(
    BOUNDARY(0, 0, -1), nx1 * ny1 * sizeof(typename MfieldsHydro::Element));
  vgrid->mp_size_send_buffer(
    BOUNDARY(0, 0, 1), nx1 * ny1 * sizeof(typename MfieldsHydro::Element));
#endif
}

// ----------------------------------------------------------------------
// vpic_define_grid

void vpic_define_grid(const Grid_t& grid)
{
#ifdef VPIC
  auto domain = grid.domain;
  auto bc = grid.bc;
  auto dt = grid.dt;

  vgrid = Grid::create();
  vgrid->setup(domain.dx, dt, grid.norm.cc, grid.norm.eps0);

  // define the grid
  define_periodic_grid(domain.corner, domain.corner + domain.length,
                       domain.gdims, domain.np);

  // set field boundary conditions
  for (int p = 0; p < grid.n_patches(); p++) {
    assert(p == 0);
    for (int d = 0; d < 3; d++) {
      bool lo = grid.atBoundaryLo(p, d);
      bool hi = grid.atBoundaryHi(p, d);

      if (lo && bc.fld_lo[d] != BND_FLD_PERIODIC) {
        Int3 bnd = {0, 0, 0};
        bnd[d] = -1;
        set_domain_field_bc(bnd, bc.fld_lo[d]);
      }

      if (hi && bc.fld_hi[d] != BND_FLD_PERIODIC) {
        Int3 bnd = {0, 0, 0};
        bnd[d] = 1;
        set_domain_field_bc(bnd, bc.fld_hi[d]);
      }
    }
  }

  // set particle boundary conditions
  for (int p = 0; p < grid.n_patches(); p++) {
    assert(p == 0);
    for (int d = 0; d < 3; d++) {
      bool lo = grid.atBoundaryLo(p, d);
      bool hi = grid.atBoundaryHi(p, d);

      if (lo && bc.prt_lo[d] != BND_PRT_PERIODIC) {
        Int3 bnd = {0, 0, 0};
        bnd[d] = -1;
        set_domain_particle_bc(bnd, bc.prt_lo[d]);
      }

      if (hi && bc.prt_hi[d] != BND_PRT_PERIODIC) {
        Int3 bnd = {0, 0, 0};
        bnd[d] = 1;
        set_domain_particle_bc(bnd, bc.prt_hi[d]);
      }
    }
  }

  grid_setup_communication();
#endif
}

// ----------------------------------------------------------------------
// vpic_define_material

#ifdef VPIC
static Material* vpic_define_material(MaterialList& material_list,
                                      const char* name, double eps,
                                      double mu = 1., double sigma = 0.,
                                      double zeta = 0.)
{
  auto m = MaterialList::create(name, eps, eps, eps, mu, mu, mu, sigma, sigma,
                                sigma, zeta, zeta, zeta);
  return material_list.append(m);
}
#else
static void vpic_define_material(MaterialList& material_list, const char* name,
                                 double eps, double mu = 1., double sigma = 0.,
                                 double zeta = 0.)
{}
#endif

// ----------------------------------------------------------------------
// vpic_define_fields

void vpic_define_fields(const Grid_t& grid)
{
#ifdef VPIC
  hydro.reset(new MfieldsHydro{grid, vgrid});
  interpolator.reset(new MfieldsInterpolator{vgrid});
  accumulator.reset(new MfieldsAccumulator{vgrid});
#endif
}

// ----------------------------------------------------------------------
// vpic_create_diagnotics

void vpic_create_diagnostics(int interval)
{
#ifdef VPIC
  diag_mixin.diagnostics_init(interval);
#endif
}

// ----------------------------------------------------------------------
// vpic_setup_diagnostics

void vpic_setup_diagnostics()
{
#ifdef VPIC
  diag_mixin.diagnostics_setup();
#endif
}

#ifdef VPIC
// ----------------------------------------------------------------------
// vpic_run_diagnostics

void vpic_run_diagnostics(Mparticles& mprts, MfieldsState& mflds)
{
  diag_mixin.diagnostics_run(mprts, mflds, *interpolator, *hydro,
                             mprts.grid().domain.np);
}
#endif
