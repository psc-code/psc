
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

// ======================================================================
// simple statistics

#define MAX_PSC_STATS 20

double psc_stats_val[MAX_PSC_STATS];
const char *psc_stats_name[MAX_PSC_STATS];
int nr_psc_stats;

int
psc_stats_register(const char *name)
{
  assert(nr_psc_stats < MAX_PSC_STATS);
  psc_stats_name[nr_psc_stats] = name;
  psc_stats_val[nr_psc_stats] = 0.;
  nr_psc_stats++;

  // we return the index + 1, so that all valid, registered tokens are > 0
  return nr_psc_stats;
}

#define psc_stats_start(n) do {				\
    psc_stats_val[n-1] = -MPI_Wtime();			\
  } while (0)

#define psc_stats_restart(n) do {			\
    psc_stats_val[n-1] -= MPI_Wtime();			\
  } while (0)

#define psc_stats_stop(n) do {				\
    psc_stats_val[n-1] += MPI_Wtime();			\
  } while (0)

#define psc_stats_val(n) psc_stats_val[n-1]

static void
psc_stats_log(struct psc *psc)
{
  int rank, size;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);

  double stats_min[nr_psc_stats], stats_max[nr_psc_stats], stats_sum[nr_psc_stats];
  MPI_Reduce(psc_stats_val, stats_min, nr_psc_stats, MPI_DOUBLE, MPI_MIN, 0, MPI_COMM_WORLD);
  MPI_Reduce(psc_stats_val, stats_max, nr_psc_stats, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
  MPI_Reduce(psc_stats_val, stats_sum, nr_psc_stats, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);

  if (rank == 0) {
    printf("    "
	   "======================================================== step %-7d ===\n",
	   psc->timestep);
    printf("    %25s %10s %10s %10s %10s\n", "name",
	   "avg", "min", "max", "total");
    for (int i = 0; i < nr_psc_stats; i++) {
      if (stats_max[i] < 1e-3)
	continue;

      printf("    %25s %10g %10g %10g %10g\n", psc_stats_name[i],
	     stats_sum[i] / size, stats_min[i], stats_max[i], stats_sum[i]);
    }
    printf("    =========================================================================\n");
  }
}

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

  static int st_nr_particles;
  static int st_time_step;
  static int st_time_particle;
  static int st_time_field;
  static int st_time_randomize;
  static int st_time_sort;
  static int st_time_collision;
  static int st_time_out_field;
  static int st_time_out_particle;
  static int st_time_balance;
  if (!st_time_step) {
    st_nr_particles = psc_stats_register("nr particles");
    st_time_step = psc_stats_register("time entire step");
    st_time_particle = psc_stats_register("time particle update");
    st_time_field = psc_stats_register("time field update");
    st_time_randomize = psc_stats_register("time randomize");
    st_time_sort = psc_stats_register("time sort");
    st_time_collision = psc_stats_register("time collision");
    st_time_out_field = psc_stats_register("time field output");
    st_time_out_particle = psc_stats_register("time particle output");
    st_time_balance = psc_stats_register("time balancing");
  }

  for (; psc->timestep < psc->prm.nmax; psc->timestep++) {
    prof_start(pr);
    psc_stats_start(st_time_step);

    psc_stats_start(st_time_out_field);
    psc_output_fields_run(psc->output_fields, psc->flds, psc->particles);
    psc_stats_stop(st_time_out_field);

    psc_stats_start(st_time_out_particle);
    psc_output_particles_run(psc->output_particles, psc->particles);
    psc_stats_stop(st_time_out_particle);

    psc_stats_start(st_time_balance);
    psc_balance_run(psc->balance, psc);
    psc_stats_stop(st_time_balance);

    psc_stats_start(st_time_randomize);
    psc_randomize_run(psc->randomize, psc->particles);
    psc_stats_stop(st_time_randomize);

    psc_stats_start(st_time_sort);
    psc_sort_run(psc->sort, psc->particles);
    psc_stats_stop(st_time_sort);

    psc_stats_start(st_time_collision);
    psc_collision_run(psc->collision, psc->particles);
    psc_stats_stop(st_time_collision);

    // field propagation n*dt -> (n+0.5)*dt
    psc_stats_start(st_time_field);
    psc_push_fields_step_a(psc->push_fields, psc->flds);
    psc_stats_stop(st_time_field);

    // particle propagation n*dt -> (n+1.0)*dt
    psc_stats_start(st_time_particle);
    psc_push_particles_run(psc->push_particles, psc->particles, psc->flds);
    // FIXME, this isn't part of particle pushing
    psc_bnd_add_ghosts(psc->bnd, psc->flds, JXI, JXI + 3);
    psc_bnd_fill_ghosts(psc->bnd, psc->flds, JXI, JXI + 3);
    psc_bnd_exchange_particles(psc->bnd, psc->particles);
    psc_stats_stop(st_time_particle);

    psc_push_photons_run(psc->mphotons);
    psc_bnd_exchange_photons(psc->bnd, psc->mphotons);
    psc_event_generator_run(psc->event_generator, psc->particles, psc->flds, psc->mphotons);

    // field propagation (n+0.5)*dt -> (n+1.0)*dt
    psc_stats_start(st_time_field);
    psc_push_fields_step_b(psc->push_fields, psc->flds);
    psc_stats_stop(st_time_field);

    // FIXME, do a mparticles func for this
    psc_stats_val(st_nr_particles) = 0;
    psc_foreach_patch(psc, p) {
      psc_stats_val(st_nr_particles) += psc->particles->p[p].n_part;
    }
    psc_stats_stop(st_time_step);

    psc_stats_log(psc);
    prof_stop(pr);
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

  psc_write_checkpoint(psc);
}
