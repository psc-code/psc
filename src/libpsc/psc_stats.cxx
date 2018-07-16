
#include "psc.h"

// ======================================================================
// psc_stats: simple statistics

double psc_stats_val[MAX_PSC_STATS + 1];
int nr_psc_stats = 1;

static const char *psc_stats_name[MAX_PSC_STATS + 1];

int
psc_stats_register(const char *name)
{
  assert(nr_psc_stats <= MAX_PSC_STATS);
  psc_stats_name[nr_psc_stats] = name;
  psc_stats_val[nr_psc_stats] = 0.;
  return nr_psc_stats++;
}

void
psc_stats_log(int timestep)
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
	   timestep);
    printf("    %25s %10s %10s %10s %10s\n", "name",
	   "avg", "min", "max", "total");
    for (int i = 1; i < nr_psc_stats; i++) {
      if (stats_max[i] < 1e-3)
	continue;

      printf("    %25s %10g %10g %10g %10g\n", psc_stats_name[i],
	     stats_sum[i] / size, stats_min[i], stats_max[i], stats_sum[i]);
    }
    printf("    =========================================================================\n");
  }

  // reset counters
  for (int i = 1; i < nr_psc_stats; i++) {
    psc_stats_val[i] = 0.;
  }
}

