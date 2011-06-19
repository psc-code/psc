
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

enum {
  STAT_NR_PARTICLES,
  STAT_TIME_STEP,
  STAT_TIME_PARTICLE,
  STAT_TIME_FIELD,
  STAT_TIME_RANDOMIZE,
  STAT_TIME_SORT,
  STAT_TIME_COLLISION,
  STAT_TIME_OUT_FIELD,
  STAT_TIME_OUT_PARTICLE,
  STAT_TIME_BALANCE,
  NR_STATS,
};

static const char *stat_name[NR_STATS] = {
  [STAT_NR_PARTICLES]      = "nr particles",
  [STAT_TIME_STEP]         = "time entire step",
  [STAT_TIME_PARTICLE]     = "time particle update",
  [STAT_TIME_FIELD]        = "time field update",
  [STAT_TIME_RANDOMIZE]    = "time part. randomize",
  [STAT_TIME_SORT]         = "time particle sort",
  [STAT_TIME_COLLISION]    = "time part. collision",
  [STAT_TIME_OUT_FIELD]    = "time field output",
  [STAT_TIME_OUT_PARTICLE] = "time particle output",
  [STAT_TIME_BALANCE]      = "time balance",
};

#define time_start(n) do {			\
    stats[n] = -MPI_Wtime();			\
  } while (0)

#define time_restart(n) do {			\
    stats[n] -= MPI_Wtime();			\
  } while (0)

#define time_stop(n) do {			\
    stats[n] += MPI_Wtime();			\
  } while (0)

static void
psc_log_step(struct psc *psc, double stats[NR_STATS])
{
  int rank, size;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);

  double stats_min[NR_STATS], stats_max[NR_STATS], stats_sum[NR_STATS];
  MPI_Reduce(stats, stats_min, NR_STATS, MPI_DOUBLE, MPI_MIN, 0, MPI_COMM_WORLD);
  MPI_Reduce(stats, stats_max, NR_STATS, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
  MPI_Reduce(stats, stats_sum, NR_STATS, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);

  if (rank == 0) {
    printf("    "
	   "======================================================== step %-7d ===\n",
	   psc->timestep);
    printf("    %25s %10s %10s %10s %10s\n", "name",
	   "avg", "min", "max", "total");
    for (int i = 0; i < NR_STATS; i++) {
      printf("    %25s %10g %10g %10g %10g\n", stat_name[i],
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

  mparticles_base_t *particles = &psc->particles;

  double stats[NR_STATS];

  for (; psc->timestep < psc->prm.nmax; psc->timestep++) {
    prof_start(pr);
    time_start(STAT_TIME_STEP);

    time_start(STAT_TIME_OUT_FIELD);
    psc_output_fields_run(psc->output_fields, psc->flds, particles);
    time_stop(STAT_TIME_OUT_FIELD);

    time_start(STAT_TIME_OUT_PARTICLE);
    psc_output_particles_run(psc->output_particles, particles);
    time_stop(STAT_TIME_OUT_PARTICLE);

    time_start(STAT_TIME_BALANCE);
    psc_balance_run(psc->balance, psc);
    time_stop(STAT_TIME_BALANCE);

    time_start(STAT_TIME_RANDOMIZE);
    psc_randomize_run(psc->randomize, particles);
    time_stop(STAT_TIME_RANDOMIZE);

    time_start(STAT_TIME_SORT);
    psc_sort_run(psc->sort, particles);
    time_stop(STAT_TIME_SORT);

    time_start(STAT_TIME_COLLISION);
    psc_collision_run(psc->collision, particles);
    time_stop(STAT_TIME_COLLISION);

    // field propagation n*dt -> (n+0.5)*dt
    time_start(STAT_TIME_FIELD);
    psc_push_fields_step_a(psc->push_fields, psc->flds);
    time_stop(STAT_TIME_FIELD);

    // particle propagation n*dt -> (n+1.0)*dt
    time_start(STAT_TIME_PARTICLE);
    psc_push_particles_run(psc->push_particles, particles, psc->flds);
    psc_bnd_add_ghosts(psc->bnd, psc->flds, JXI, JXI + 3);
    psc_bnd_fill_ghosts(psc->bnd, psc->flds, JXI, JXI + 3);
    psc_bnd_exchange_particles(psc->bnd, particles);
    time_stop(STAT_TIME_PARTICLE);

    psc_push_photons_run(&psc->mphotons);
    psc_bnd_exchange_photons(psc->bnd, &psc->mphotons);
    psc_event_generator_run(psc->event_generator, particles, psc->flds, &psc->mphotons);

    // field propagation (n+0.5)*dt -> (n+1.0)*dt
    time_restart(STAT_TIME_FIELD);
    psc_push_fields_step_b(psc->push_fields, psc->flds);
    time_stop(STAT_TIME_FIELD);

    stats[STAT_NR_PARTICLES] = 0;
    psc_foreach_patch(psc, p) {
      stats[STAT_NR_PARTICLES] += particles->p[p].n_part;
    }
    time_stop(STAT_TIME_STEP);
    psc_log_step(psc, stats);
    // FIXME, check whether cpu time expired?
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
