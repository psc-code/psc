
#pragma once

// ----------------------------------------------------------------------
// psc_stats: simple statistics

#define MAX_PSC_STATS 20

extern double psc_stats_val[MAX_PSC_STATS+1]; // [0] is left empty
extern int nr_psc_stats;

int psc_stats_register(const char *name);
void psc_stats_log(int timestep);

#define psc_stats_start(n) do {				\
    psc_stats_val[n] -= MPI_Wtime();			\
  } while (0)

#define psc_stats_stop(n) do {				\
    psc_stats_val[n] += MPI_Wtime();			\
  } while (0)

// These are general statistics categories to be used in different parts
// of the code as appropriate.

extern int st_time_output;   //< time spent in output
extern int st_time_comm;     //< time spent in communications
extern int st_time_particle; //< time spent in particle computation
extern int st_time_field;    //< time spent in field computation

