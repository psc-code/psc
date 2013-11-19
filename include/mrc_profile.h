
#ifndef PROFILE_H
#define PROFILE_H

#include <assert.h>
#include <stdio.h>
#include <mpi.h>

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

//#define PROF_UOPS
#define PROF_CACHE

#ifdef HAVE_LIBPAPI

#if defined(PROF_UOPS) || defined(PROF_CACHE)
#define NR_EVENTS (4)
#else
#define NR_EVENTS (3)
#endif

#else
#define NR_EVENTS (0)
#endif

struct prof_info {
  int cnt;
  long long time;
  long long counters[NR_EVENTS];
};

#define MAX_PROF (100)

extern struct prof_globals {
  int event_set;
  struct prof_info info[MAX_PROF];
} prof_globals;

#ifdef HAVE_LIBPAPI

#include <papi.h>
#include <stdlib.h>
static inline void
prof_start(int pr)
{
  int rc;
  long long counters[NR_EVENTS];

  pr--;
  assert(pr < MAX_PROF);
  prof_globals.info[pr].time -= PAPI_get_real_usec();
  rc = PAPI_read(prof_globals.event_set, counters);
  assert(rc == PAPI_OK);

  for (int i = 0; i < NR_EVENTS; i++) {
    prof_globals.info[pr].counters[i] -= counters[i];
  }
}

static inline void
prof_stop(int pr)
{
  int rc;
  long long counters[NR_EVENTS];
  long long t1 = PAPI_get_real_usec();

  pr--;
  rc = PAPI_read(prof_globals.event_set, counters);
  assert(rc == PAPI_OK);

  for (int i = 0; i < NR_EVENTS; i++) {
    prof_globals.info[pr].counters[i] += counters[i];
  }
  
  prof_globals.info[pr].time += t1;
  prof_globals.info[pr].cnt++;
}

#else // ! HAVE_LIBPAPI

#include <sys/time.h>
#include <stdlib.h>

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

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

#endif

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
