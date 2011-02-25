
#include "mrc_profile.h"

#include <assert.h>
#include <string.h>

struct prof_data {
  const char *name;
  float simd;
  int flops;
  int bytes;
};

static int prof_inited;
static int nr_prof_data;
static struct prof_data prof_data[MAX_PROF];

#ifdef HAVE_LIBPAPI

#include <papi.h>

static const int events[] = {
  PAPI_TOT_CYC, // this one should be kept unchanged
  PAPI_TOT_INS, 
#ifdef PROF_UOPS
  0x4000005b,
  0x4000205b,
#elif defined(PROF_CACHE)
  //  0x40040037, // prefetch_hit
  //  0x40080037, // prefetch_miss
  //  0x40008037, // ld_hit
  //  0x40010037, // ld_miss
  0x40000006, // llc_misses
  //  0x4000082f, // l1d_prefetch:requests
  //  0x4000102f, // l1d_prefetch:triggers
  //  0x40000417, // dtlb_misses
  PAPI_FP_OPS,
#else
  PAPI_FP_OPS,
#endif
};

void
prof_init(void)
{
  int rc;

  if (prof_inited)
    return;

  // assert(ARRAYSIZE(events) == NR_EVENTS);

  rc = PAPI_library_init(PAPI_VER_CURRENT);
  assert(rc == PAPI_VER_CURRENT);

  memset(&prof_globals, 0, sizeof(prof_globals));
  prof_globals.event_set = PAPI_NULL;
  rc = PAPI_create_eventset(&prof_globals.event_set);
  assert(rc == PAPI_OK);
  
  for (int i = 0; i < NR_EVENTS; i++) {
    rc = PAPI_add_event(prof_globals.event_set, events[i]);
    assert(rc == PAPI_OK);
  }

  rc = PAPI_start(prof_globals.event_set);
  assert(rc == PAPI_OK);

  prof_inited = 1;
}

void
prof_print_file(FILE *f)
{
  int rc;

  fprintf(f, "%20s %7s %4s %7s", "", "tottime", "cnt", "time");
  for (int i = 0; i < NR_EVENTS; i++) {
    char name[200];
    rc = PAPI_event_code_to_name(events[i], name);
    assert(rc == PAPI_OK);
    fprintf(f, " %12s", name);
  }
  fprintf(f, " %12s", "FLOPS");
  fprintf(f, " %12s", "MBytes");
  fprintf(f, "\n");

  for (int pr = 0; pr < nr_prof_data; pr++) {
    struct prof_info *pinfo = &prof_globals.info[pr];
    int cnt = pinfo->cnt;
    if (!cnt)
      continue;

    long long cycles = pinfo->counters[0];
    float rtime = pinfo->time;
    fprintf(f, "%20s %7g %4d %7g", prof_data[pr].name, rtime/1e3, cnt, 
	    rtime / 1e3 / cnt);
    for (int i = 0; i < NR_EVENTS; i++) {
      long long counter = pinfo->counters[i];
      if (events[i] == PAPI_FP_OPS) {
	counter *= prof_data[pr].simd;
      }
      fprintf(f, " %12lld", counter / cnt);
    }
    fprintf(f, " %12d", prof_data[pr].flops);
    fprintf(f, " %12g", prof_data[pr].bytes / 1e6);
    fprintf(f, "\n");

    fprintf(f, "%20s %7s %4s %7s", "", "", "", "");
    for (int i = 0; i < NR_EVENTS; i++) {
      long long counter = pinfo->counters[i];
      if (events[i] == PAPI_FP_OPS) {
	counter *= prof_data[pr].simd;
      }
      fprintf(f, " %12g", (float) counter / cycles);
    }
    fprintf(f, " %12g", (float) prof_data[pr].flops / (cycles/cnt));
    fprintf(f, " %12g", (float) prof_data[pr].bytes / (cycles/cnt));
    fprintf(f, " / cycle\n");

    fprintf(f, "%20s %7s %4s %7s", "", "", "", "");
    for (int i = 0; i < NR_EVENTS; i++) {
      long long counter = pinfo->counters[i];
      if (events[i] == PAPI_FP_OPS) {
	counter *= prof_data[pr].simd;
      }
      fprintf(f, " %12g", counter / rtime);
    }
    fprintf(f, " %12g", (float) prof_data[pr].flops / (rtime/cnt));
    fprintf(f, " %12g", (float) prof_data[pr].bytes / (rtime/cnt));
    fprintf(f, " / usec\n");
  }
}

#else // !LIBPAPI

void
prof_init(void)
{
  prof_inited = 1;
}

void
prof_print_file(FILE *f)
{
  fprintf(f, "%19s %7s %4s %7s", "", "tottime", "cnt", "time");
  fprintf(f, " %12s", "FLOPS");
  fprintf(f, " %12s", "MFLOPS/sec");
  fprintf(f, " %12s", "MBytes/sec");
  fprintf(f, "\n");

  for (int pr = 0; pr < nr_prof_data; pr++) {
    struct prof_info *pinfo = &prof_globals.info[pr];
    int cnt = pinfo->cnt;
    float rtime = pinfo->time;
    if (!cnt || rtime == 0.)
      continue;

    fprintf(f, "%-19s %7g %4d %7g", prof_data[pr].name, rtime/1e3, cnt, 
	   rtime / 1e3 / cnt);
    fprintf(f, " %12d", prof_data[pr].flops);
    fprintf(f, " %12g", (float) prof_data[pr].flops / (rtime/cnt));
    fprintf(f, " %12g", prof_data[pr].bytes / (rtime/cnt));
    fprintf(f, "\n");
  }
}

#endif

int
prof_register(const char *name, float simd, int flops, int bytes)
{
  if (!prof_inited) {
    prof_init();
  }

  assert(nr_prof_data < MAX_PROF);
  struct prof_data *p = &prof_data[nr_prof_data++];

  p->name = name;
  p->simd = simd;
  p->flops = flops;
  p->bytes = bytes;

  return nr_prof_data;
}

void
prof_print()
{
  prof_print_file(stdout);
}

void
prof_print_mpi(MPI_Comm comm)
{
  float times[nr_prof_data], times_avg[nr_prof_data];
  float times_min[nr_prof_data], times_max[nr_prof_data];

  for (int pr = 0; pr < nr_prof_data; pr++) {
    struct prof_info *pinfo = &prof_globals.info[pr];
    if (pinfo->cnt > 0) {
      times[pr] = pinfo->time / pinfo->cnt / 1e3;
    } else {
      times[pr] = -1;
    }
  }

  int rank, size;
  MPI_Comm_rank(comm, &rank);
  MPI_Comm_size(comm, &size);

  MPI_Reduce(times, times_avg, nr_prof_data, MPI_FLOAT, MPI_SUM, 0, comm);
  MPI_Reduce(times, times_min, nr_prof_data, MPI_FLOAT, MPI_MIN, 0, comm);
  MPI_Reduce(times, times_max, nr_prof_data, MPI_FLOAT, MPI_MAX, 0, comm);

  if (rank == 0) {
    printf("%19s %10s %10s %10s\n", "", "avg", "min", "max");
    for (int pr = 0; pr < nr_prof_data; pr++) {
      times_avg[pr] /= size;
      
      printf("%-19s %10.2f %10.2f %10.2f\n", prof_data[pr].name, times_avg[pr],
	     times_min[pr], times_max[pr]);
    }
  }

  for (int pr = 0; pr < nr_prof_data; pr++) {
    struct prof_info *pinfo = &prof_globals.info[pr];
    pinfo->time = 0.;
    pinfo->cnt = 0;
  }
}
