
#include "psc.h"

// ======================================================================
// simple statistics

enum {
  STAT_NR_PARTICLES,
  STAT_TIME_STEP,
  STAT_TIME_PARTICLE,
  STAT_TIME_FIELD,
  STAT_TIME_OUT_FIELD,
  STAT_TIME_OUT_PARTICLE,
  NR_STATS,
};

static const char *stat_name[NR_STATS] = {
  [STAT_NR_PARTICLES]      = "nr particles",
  [STAT_TIME_STEP]         = "time entire step",
  [STAT_TIME_PARTICLE]     = "time particle update",
  [STAT_TIME_FIELD]        = "time field update",
  [STAT_TIME_OUT_FIELD]    = "time field output",
  [STAT_TIME_OUT_PARTICLE] = "time particle output",
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
psc_log_step(double stats[NR_STATS])
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
	   psc.timestep);
    printf("    %25s %10s %10s %10s %10s\n", "name",
	   "avg", "min", "max", "total");
    for (int i = 0; i < NR_STATS; i++) {
      printf("    %25s %10g %10g %10g %10g\n", stat_name[i],
	     stats_sum[i] / size, stats_min[i], stats_max[i], stats_sum[i]);
    }
    printf("    =========================================================================\n");
  }
}

#define ALLOC_field_fortran_F77 F77_FUNC_(alloc_field_fortran, ALLOC_FIELD_FORTRAN)
#define SETUP_field_F77 F77_FUNC_(setup_field, SETUP_FIELD)

void ALLOC_field_fortran_F77(void);
void SETUP_field_F77(void);

/////////////////////////////////////////////////////////////////////////
/// Main time integration loop.
///

void
psc_integrate()
{
  ALLOC_field_fortran_F77();
  SETUP_field_F77();

  double stats[NR_STATS];

  for (; psc.timestep < psc.prm.nmax; psc.timestep++) {
    time_start(STAT_TIME_STEP);

    time_start(STAT_TIME_OUT_FIELD);
    psc_out_field();
    time_stop(STAT_TIME_OUT_FIELD);

    time_start(STAT_TIME_OUT_PARTICLE);
    psc_out_particles();
    time_stop(STAT_TIME_OUT_PARTICLE);

    // field propagation n*dt -> (n+0.5)*dt
    time_start(STAT_TIME_FIELD);
    psc_push_field_a();
    time_stop(STAT_TIME_FIELD);

    // particle propagation n*dt -> (n+1.0)*dt
    time_start(STAT_TIME_PARTICLE);
    psc_push_particles();
    int flds[] = { JXI, JYI, JZI, -1 };
    for (int i = 0; flds[i] >= 0; i++) {
      int m = flds[i];
      psc_add_ghosts(m);
      psc_fill_ghosts(m);
    }
    psc_exchange_particles();
    time_stop(STAT_TIME_PARTICLE);

    // field propagation (n+0.5)*dt -> (n+1.0)*dt
    time_restart(STAT_TIME_FIELD);
    psc_push_field_b();
    time_stop(STAT_TIME_FIELD);

    stats[STAT_NR_PARTICLES] = psc.n_part;
    time_stop(STAT_TIME_STEP);
    psc_log_step(stats);
    // FIXME, check whether cpu time expired?
  }
}
