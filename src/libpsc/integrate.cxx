
#include "psc.h"
#include "psc_method.h"
#include "psc_push_particles.h"
#include "psc_push_fields.h"
#include "psc_marder.h"
#include "psc_bnd_particles.h"
#include "psc_collision.h"
#include "psc_randomize.h"
#include "psc_sort.h"
#include "psc_event_generator.h"
#include "psc_balance.h"
#include "psc_checks.h"
#include "particles.hxx"
#include "push_particles.hxx"
#include "push_fields.hxx"

#include <mrc_common.h>
#include <mrc_profile.h>

int st_time_output;
int st_time_comm;
int st_time_particle;
int st_time_field;

#define psc_ops(psc) ((struct psc_ops *)((psc)->obj.ops))

// ----------------------------------------------------------------------
// psc_print_profiling

void
psc_print_profiling(struct psc *psc)
{
  int size;
  MPI_Comm_size(MPI_COMM_WORLD, &size);
  if (size > 1 && !psc->prm.detailed_profiling) {
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
// psc_output

void
psc_output(struct psc *psc)
{
  if (psc_ops(psc) && psc_ops(psc)->output) {
    return psc_ops(psc)->output(psc);
  }

  psc_method_output(psc->method, psc);
}

/////////////////////////////////////////////////////////////////////////
/// psc_step
///

// This measures the time spent pushing particles and fields, exclusive of
// communication.
// Only works correctly for push_fields "variant 1"!
int pr_time_step_no_comm; // FIXME, don't like globals

void
psc_step(struct psc *psc)
{
  if (psc_ops(psc) && psc_ops(psc)->step) {
    psc_ops(psc)->step(psc);
    return;
  }

  if (!pr_time_step_no_comm) {
    pr_time_step_no_comm = prof_register("time step w/o comm", 1., 0, 0);
  }

  // default psc_step() implementation

  mparticles_base_t mprts(psc->particles);
  mfields_base_t mflds(psc->flds);
  PscPushParticlesBase<void> pushp(psc->push_particles);
  PscPushFields<void> pushf(psc->push_fields);
#if 0
  mpi_printf(psc_comm(psc), "**** Step %d / %d, Time %g\n", psc->timestep + 1,
	     psc->prm.nmax, psc->timestep * psc->dt);
#endif

  // x^{n+1/2}, p^{n}, E^{n+1/2}, B^{n+1/2}

  psc_balance_run(psc->balance, psc);

  psc_randomize_run(psc->randomize, psc->particles);
  psc_sort_run(psc->sort, psc->particles);

  prof_start(pr_time_step_no_comm);
  prof_stop(pr_time_step_no_comm); // actual measurements are done w/ restart

  psc_collision_run(psc->collision, psc->particles);
  
  //psc_bnd_particles_open_calc_moments(psc->bnd_particles, psc->particles);

  psc_checks_continuity_before_particle_push(psc->checks, psc);

  // particle propagation p^{n} -> p^{n+1}, x^{n+1/2} -> x^{n+3/2}
  pushp(mprts, mflds);
  // x^{n+3/2}, p^{n+1}, E^{n+1/2}, B^{n+1/2}, j^{n+1}
    
  // field propagation B^{n+1/2} -> B^{n+1}
  pushf.advance_H(mflds, .5);
  // x^{n+3/2}, p^{n+1}, E^{n+1/2}, B^{n+1}, j^{n+1}

  psc_bnd_particles_exchange(psc->bnd_particles, psc->particles);
  
  psc_event_generator_run(psc->event_generator, psc->particles, psc->flds);
  
  // field propagation E^{n+1/2} -> E^{n+3/2}
  pushf.advance_b2(mflds);
  // x^{n+3/2}, p^{n+1}, E^{n+3/2}, B^{n+1}

  // field propagation B^{n+1} -> B^{n+3/2}
  pushf.advance_a(mflds);
  // x^{n+3/2}, p^{n+1}, E^{n+3/2}, B^{n+3/2}

  psc_checks_continuity_after_particle_push(psc->checks, psc);

  // E at t^{n+3/2}, particles at t^{n+3/2}
  // B at t^{n+3/2} (Note: that is not it's natural time,
  // but div B should be == 0 at any time...)
  psc_marder_run(psc->marder, psc->flds, psc->particles);
    
  psc_checks_gauss(psc->checks, psc);

  psc_push_particles_prep(psc->push_particles, psc->particles, psc->flds);
}

/////////////////////////////////////////////////////////////////////////
/// Main time integration loop.
///

void
psc_integrate(struct psc *psc)
{
  psc_method_initialize(psc->method, psc);
  mpi_printf(psc_comm(psc), "Initialization complete.\n");
  
  static int pr;
  if (!pr) {
    pr = prof_register("psc_step", 1., 0, 0);
  }

  int st_nr_particles = psc_stats_register("nr particles");
  int st_time_step = psc_stats_register("time entire step");

  // generic stats categories
  st_time_particle = psc_stats_register("time particle update");
  st_time_field = psc_stats_register("time field update");
  st_time_comm = psc_stats_register("time communication");
  st_time_output = psc_stats_register("time output");

  mpi_printf(psc_comm(psc), "*** Advancing\n");
  double elapsed = MPI_Wtime();

  bool first_iteration = true;
  while (psc->timestep < psc->prm.nmax) {
    prof_start(pr);
    psc_stats_start(st_time_step);

    if (!first_iteration &&
	psc->prm.write_checkpoint_every_step > 0 &&
	psc->timestep % psc->prm.write_checkpoint_every_step == 0) {
      psc_write_checkpoint(psc);
    }
    first_iteration = false;

    psc_step(psc);

    psc->timestep++; // FIXME, too hacky
    psc_output(psc);
    
    psc_stats_stop(st_time_step);
    prof_stop(pr);

    mparticles_base_t mprts(psc->particles);
    psc_stats_val[st_nr_particles] = mprts->get_n_prts();

    if (psc->timestep % psc->prm.stats_every == 0) {
      psc_stats_log(psc);
      psc_print_profiling(psc);
    }

    if (psc->prm.wallclock_limit > 0.) {
      double wallclock_elapsed = MPI_Wtime() - psc->time_start;
      double wallclock_elapsed_max;
      MPI_Allreduce(&wallclock_elapsed, &wallclock_elapsed_max, 1, MPI_DOUBLE, MPI_MAX,
		    MPI_COMM_WORLD);
      
      if (wallclock_elapsed_max > psc->prm.wallclock_limit) {
	mpi_printf(MPI_COMM_WORLD, "WARNING: Max wallclock time elapsed!\n");
	break;
      }
    }
  }

  if (psc->prm.write_checkpoint) {
    psc_write_checkpoint(psc);
  }

  // FIXME, merge with existing handling of wallclock time
  elapsed = MPI_Wtime() - elapsed;

  int  s = (int)elapsed, m  = s/60, h  = m/60, d  = h/24, w = d/ 7;
  /**/ s -= m*60,        m -= h*60, h -= d*24, d -= w*7;
  mpi_printf(psc_comm(psc), "*** Finished (%gs / %iw:%id:%ih:%im:%is elapsed)\n",
	     elapsed, w, d, h, m, s );
  
}

