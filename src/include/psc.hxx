
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

struct MaterialList;

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

    for (int d = 0; d < 3; d++) {
      if (grid.isInvar(d) != Dim::is_invar(d)) {
        LOG_ERROR("dimension %d is%s invariant, but gdims[%d]=%d\n", d,
                  Dim::is_invar(d) ? "" : " not", d, grid.domain.gdims[d]);
      }
    }

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
    bndf_.fill_ghosts_H(mflds_);
    bnd_.fill_ghosts(mflds_, HX, HX + 3);

    bnd_.fill_ghosts(mflds_, JXI, JXI + 3);

    bndf_.fill_ghosts_E(mflds_);
    bnd_.fill_ghosts(mflds_, EX, EX + 3);

    // FIXME: do a half-step on p to bring it to its natural time,
    // p^{n+1/2} -> p^{n+1}
    // pushp_.stagger(mprts, mflds);

    if (checks_.gauss.should_do_check(0)) {
      mpi_printf(grid().comm(), "Checking gauss.\n");
      checks_.gauss(mprts_, mflds_);
    }

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

  // ----------------------------------------------------------------------
  // step

  void step()
  {
    using Dim = typename PscConfig::Dim;

    static int pr_sort, pr_collision, pr_checks, pr_push_prts, pr_push_flds,
      pr_bndp, pr_bndf, pr_marder, pr_inject_prts, pr_external_current;
    if (!pr_sort) {
      pr_sort = prof_register("step_sort", 1., 0, 0);
      pr_collision = prof_register("step_collision", 1., 0, 0);
      pr_push_prts = prof_register("step_push_prts", 1., 0, 0);
      pr_inject_prts = prof_register("step_inject_prts", 1., 0, 0);
      pr_external_current = prof_register("step_external_current", 1., 0, 0);
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

    if (checks_.continuity.should_do_check(timestep)) {
      mpi_printf(comm, "***** Checking continuity (1 of 2)...\n");
      prof_start(pr_checks);
      checks_.continuity.before_particle_push(mprts_);
      prof_stop(pr_checks);
    }

    // === particle propagation p^{n} -> p^{n+1}, x^{n+1/2} -> x^{n+3/2}
    mpi_printf(comm, "***** Push particles...\n");
    prof_start(pr_push_prts);
    pushp_.push_mprts(mprts_, mflds_);
    prof_stop(pr_push_prts);

    // === particle injection
    prof_start(pr_inject_prts);
    inject_particles_(mprts_, mflds_);
    prof_stop(pr_inject_prts);

    // === external current
    prof_start(pr_external_current);
    this->ext_current_(grid(), mflds_);
    prof_stop(pr_external_current);

    mpi_printf(comm, "***** Bnd particles...\n");
    prof_start(pr_bndp);
    bndp_(mprts_);
    prof_stop(pr_bndp);

    mpi_printf(comm, "***** Bnd fields J...\n");
    prof_start(pr_bndf);
    bndf_.add_ghosts_J(mflds_);
    bnd_.add_ghosts(mflds_, JXI, JXI + 3);
    bnd_.fill_ghosts(mflds_, JXI, JXI + 3);
    prof_stop(pr_bndf);
    // state is now: x^{n+3/2}, p^{n+1}, E^{n+1/2}, B^{n+1/2}, j^{n+1}

    // === field propagation B^{n+1/2} -> B^{n+1}
    mpi_printf(comm, "***** Push fields B (1 of 2)...\n");
    prof_start(pr_push_flds);
    pushf_.push_H(mflds_, .5, Dim{});
    prof_stop(pr_push_flds);

    mpi_printf(comm, "***** Bnd fields B (1 of 2)...\n");
    prof_start(pr_bndf);
    bndf_.fill_ghosts_H(mflds_);
    bnd_.fill_ghosts(mflds_, HX, HX + 3);
    prof_stop(pr_bndf);
    // state is now: x^{n+3/2}, p^{n+1}, E^{n+1/2}, B^{n+1}, j^{n+1}

    // === field propagation E^{n+1/2} -> E^{n+3/2}
    mpi_printf(comm, "***** Push fields E...\n");
    prof_restart(pr_push_flds);
    pushf_.push_E(mflds_, 1., Dim{});
    prof_stop(pr_push_flds);

    mpi_printf(comm, "***** Bnd fields E...\n");
    prof_restart(pr_bndf);
    bndf_.fill_ghosts_E(mflds_);
    bnd_.fill_ghosts(mflds_, EX, EX + 3);
    prof_stop(pr_bndf);
    // state is now: x^{n+3/2}, p^{n+1}, E^{n+3/2}, B^{n+1}

    // === field propagation B^{n+1} -> B^{n+3/2}
    mpi_printf(comm, "***** Push fields B (2 of 2)...\n");
    prof_restart(pr_push_flds);
    pushf_.push_H(mflds_, .5, Dim{});
    prof_stop(pr_push_flds);

    mpi_printf(comm, "***** Bnd fields B (2 of 2)...\n");
    prof_start(pr_bndf);
    bndf_.fill_ghosts_H(mflds_);
    bnd_.fill_ghosts(mflds_, HX, HX + 3);
    prof_stop(pr_bndf);
    // state is now: x^{n+3/2}, p^{n+1}, E^{n+3/2}, B^{n+3/2}

    if (checks_.continuity.should_do_check(timestep)) {
      mpi_printf(comm, "***** Checking continuity (2 of 2)...\n");
      prof_restart(pr_checks);
      checks_.continuity.after_particle_push(mprts_, mflds_);
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

    if (checks_.gauss.should_do_check(timestep)) {
      mpi_printf(comm, "***** Checking gauss...\n");
      prof_restart(pr_checks);
      checks_.gauss(mprts_, mflds_);
      prof_stop(pr_checks);
    }

    // psc_push_particles_prep(psc->push_particles, psc->particles, psc->flds);
  }

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

  // ----------------------------------------------------------------------
  // diagnostics

  void diagnostics() { diagnostics_(mprts_, mflds_); }

  // ----------------------------------------------------------------------
  // print_status

  void print_status()
  {
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
  template <typename Mparticles, typename MfieldsState>
  void operator()(Mparticles& mprts, MfieldsState& mflds)
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
{}

// ----------------------------------------------------------------------
// set_domain_field_bc

void set_domain_field_bc(Int3 bnd, int bc) {}

// ----------------------------------------------------------------------
// set_domain_particle_bc

void set_domain_particle_bc(Int3 bnd, int bc) {}

void grid_setup_communication() {}

// ----------------------------------------------------------------------
// vpic_define_grid

void vpic_define_grid(const Grid_t& grid) {}

// ----------------------------------------------------------------------
// vpic_define_material

static void vpic_define_material(MaterialList& material_list, const char* name,
                                 double eps, double mu = 1., double sigma = 0.,
                                 double zeta = 0.)
{}

// ----------------------------------------------------------------------
// vpic_define_fields

void vpic_define_fields(const Grid_t& grid) {}

// ----------------------------------------------------------------------
// vpic_create_diagnotics

void vpic_create_diagnostics(int interval) {}

// ----------------------------------------------------------------------
// vpic_setup_diagnostics

void vpic_setup_diagnostics() {}
