
#ifndef PROFILE_H
#define PROFILE_H

#include <mrc_config.h>

#include <assert.h>
#include <mpi.h>
#include <stdio.h>

#ifdef HAVE_NVTX
#include <nvToolsExt.h>
#endif

#define NR_EVENTS (0)

struct prof_info
{
  int cnt;
  long long time;
  long long counters[NR_EVENTS];
  int total_cnt;
  long long total_time;
};

struct prof_data
{
  const char* name;
  float simd;
  int flops;
  int bytes;
};

#define MAX_PROF (100)

extern struct prof_data prof_data[MAX_PROF];

extern struct prof_globals
{
  int event_set;
  struct prof_info info[MAX_PROF];
} prof_globals;

#include <stdlib.h>
#include <sys/time.h>

static inline void prof_start(int pr)
{
  pr--;
  assert(pr < MAX_PROF);

  struct timeval tv;
  gettimeofday(&tv, NULL);
  prof_globals.info[pr].time -= tv.tv_sec * 1000000ll + tv.tv_usec;
#ifdef HAVE_NVTX
  nvtxRangePush(prof_data[pr].name);
#endif
}

static inline void prof_restart(int pr)
{
  pr--;
  assert(pr < MAX_PROF);

  struct timeval tv;
  gettimeofday(&tv, NULL);
  prof_globals.info[pr].time -= tv.tv_sec * 1000000ll + tv.tv_usec;
  prof_globals.info[pr].cnt--;
#ifdef HAVE_NVTX
  nvtxRangePush(prof_data[pr].name);
#endif
}

static inline void prof_stop(int pr)
{
  pr--;
  assert(pr < MAX_PROF);

  struct timeval tv;
  gettimeofday(&tv, NULL);
  prof_globals.info[pr].time += tv.tv_sec * 1000000ll + tv.tv_usec;
  prof_globals.info[pr].cnt++;
#ifdef HAVE_NVTX
  nvtxRangePop();
#endif
}

#ifdef __cplusplus
#define EXTERN_C extern "C"
#else
#define EXTERN_C
#endif

EXTERN_C void prof_init(void);
EXTERN_C int prof_register(const char* name, float simd, int flops, int bytes);
EXTERN_C void prof_print(void);
EXTERN_C void prof_print_file(FILE* f);
EXTERN_C void prof_print_mpi(MPI_Comm comm);

#define prof_barrier(str)                                                      \
  do {                                                                         \
    static int pr;                                                             \
    if (!pr) {                                                                 \
      pr = prof_register(str, 1., 0, 0);                                       \
    }                                                                          \
    prof_start(pr);                                                            \
    MPI_Barrier(MPI_COMM_WORLD);                                               \
    prof_stop(pr);                                                             \
  } while (0)

#endif
