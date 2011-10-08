
#include "psc.h"
#include "psc_push_particles.h"
#include "psc_push_fields.h"
#include "psc_bnd.h"
#include "psc_collision.h"
#include "psc_randomize.h"
#include "psc_sort.h"
#include "psc_output_fields.h"
#include "psc_output_particles.h"
#include "psc_event_generator.h"
#include "psc_balance.h"

#include <mrc_common.h>
#include <mrc_profile.h>

int st_time_output;
int st_time_comm;
int st_time_particle;
int st_time_field;

#define psc_ops(psc) ((struct psc_ops *)((psc)->obj.ops))

/////////////////////////////////////////////////////////////////////////
/// psc_output
///

void
psc_output(struct psc *psc)
{
  if (psc_ops(psc) && psc_ops(psc)->output) {
    psc_ops(psc)->output(psc);
    return;
  }

  // default psc_output() implementation

  psc_output_fields_run(psc->output_fields, psc->flds, psc->particles);
  psc_output_particles_run(psc->output_particles, psc->particles);
}

/////////////////////////////////////////////////////////////////////////
/// psc_step
///

void
psc_step(struct psc *psc)
{
  if (psc_ops(psc) && psc_ops(psc)->step) {
    psc_ops(psc)->step(psc);
    return;
  }

  // default psc_step() implementation

  if (psc->use_dynamic_patches) {
    psc_patchmanager_timestep(&psc->patchmanager);
  }

  psc_output(psc);
  psc_balance_run(psc->balance, psc);
  
  psc_randomize_run(psc->randomize, psc->particles);
  psc_sort_run(psc->sort, psc->particles);
  psc_collision_run(psc->collision, psc->particles);
  
  // field propagation n*dt -> (n+0.5)*dt
  psc_push_fields_step_a(psc->push_fields, psc->flds);
  
  // particle propagation n*dt -> (n+1.0)*dt
  psc_push_particles_run(psc->push_particles, psc->particles, psc->flds);
  psc_bnd_exchange_particles(psc->bnd, psc->particles);
  
  psc_push_photons_run(psc->mphotons);
  psc_bnd_exchange_photons(psc->bnd, psc->mphotons);
  psc_event_generator_run(psc->event_generator, psc->particles, psc->flds, psc->mphotons);
  
  psc_bnd_add_ghosts(psc->bnd, psc->flds, JXI, JXI + 3);
  psc_bnd_fill_ghosts(psc->bnd, psc->flds, JXI, JXI + 3);
  
  // field propagation (n+0.5)*dt -> (n+1.0)*dt
  psc_push_fields_step_b(psc->push_fields, psc->flds);
}

extern void dynamicwindow_timestep();

/////////////////////////////////////////////////////////////////////////
/// Main time integration loop.
///

void
psc_integrate(struct psc *psc)
{
  static int pr;
  if (!pr) {
    pr = prof_register("psc_step", 1., 0, 0);
  }

  int st_nr_particles = psc_stats_register("nr particles");
  int st_nr_photons = psc_stats_register("nr photons");
  int st_time_step = psc_stats_register("time entire step");

  // generic stats categories
  st_time_particle = psc_stats_register("time particle update");
  st_time_field = psc_stats_register("time field update");
  st_time_comm = psc_stats_register("time communication");
  st_time_output = psc_stats_register("time output");

  for (; psc->timestep < psc->prm.nmax; psc->timestep++) {
    prof_start(pr);
    psc_stats_start(st_time_step);

    psc_step(psc);

    psc_stats_stop(st_time_step);
    prof_stop(pr);

    // FIXME, do a mparticles func for this
    psc_foreach_patch(psc, p) {
      psc_stats_val[st_nr_particles] += psc_mparticles_base_nr_particles_by_patch(psc->particles, p);
      psc_stats_val[st_nr_photons] += psc->mphotons->p[p].nr;
    }

    psc_stats_log(psc);
    prof_print_mpi(MPI_COMM_WORLD);

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
}
