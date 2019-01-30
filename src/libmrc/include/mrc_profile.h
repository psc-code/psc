
#ifndef PROFILE_H
#define PROFILE_H

#include <assert.h>
#include <stdio.h>
#include <mpi.h>

#include <mrc_config.h>

#define NR_EVENTS (0)

struct prof_info {
  int cnt;
  long long time;
  long long counters[NR_EVENTS];
  int total_cnt;
  long long total_time;
};

#define MAX_PROF (100)

extern struct prof_globals {
  int event_set;
  struct prof_info info[MAX_PROF];
} prof_globals;

#include <sys/time.h>
#include <stdlib.h>

static inline void
prof_start(int pr)
{
  pr--;
  assert(pr < MAX_PROF);

  struct timeval tv;
  gettimeofday(&tv, NULL);
  prof_globals.info[pr].time -= tv.tv_sec * 1000000ll + tv.tv_usec;
}

static inline void
prof_restart(int pr)
{
  pr--;
  assert(pr < MAX_PROF);

  struct timeval tv;
  gettimeofday(&tv, NULL);
  prof_globals.info[pr].time -= tv.tv_sec * 1000000ll + tv.tv_usec;
  prof_globals.info[pr].cnt--;
}

static inline void
prof_stop(int pr)
{
  pr--;
  assert(pr < MAX_PROF);

  struct timeval tv;
  gettimeofday(&tv, NULL);
  prof_globals.info[pr].time += tv.tv_sec * 1000000ll + tv.tv_usec;
  prof_globals.info[pr].cnt++;
}

#ifdef __cplusplus
#define EXTERN_C extern "C"
#else
#define EXTERN_C
#endif

EXTERN_C void prof_init(void);
EXTERN_C int  prof_register(const char *name, float simd, int flops, int bytes);
EXTERN_C void prof_print(void);
EXTERN_C void prof_print_file(FILE *f);
EXTERN_C void prof_print_mpi(MPI_Comm comm);

#endif
