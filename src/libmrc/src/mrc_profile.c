
#include "mrc_profile.h"

#include <assert.h>
#include <string.h>

struct prof_globals prof_globals;

struct prof_data {
  const char *name;
  float simd;
  int flops;
  int bytes;
};

static int prof_inited;
static int nr_prof_data;
static struct prof_data prof_data[MAX_PROF];

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
  for (int pr = 0; pr < nr_prof_data; pr++) {
    struct prof_info *pinfo = &prof_globals.info[pr];
    pinfo->time = 0.;
    pinfo->cnt = 0;
  }
}

void
prof_print_mpi(MPI_Comm comm)
{
  float times[nr_prof_data], times_avg[nr_prof_data];
  float times_min[nr_prof_data], times_max[nr_prof_data];

  for (int pr = 0; pr < nr_prof_data; pr++) {
    struct prof_info *pinfo = &prof_globals.info[pr];
    if (pinfo->cnt > 0) {
      times[pr] = pinfo->time;
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

  for (int pr = 0; pr < nr_prof_data; pr++) {
    struct prof_info *pinfo = &prof_globals.info[pr];
    pinfo->total_time += pinfo->time;
    pinfo->total_cnt += pinfo->cnt;
  }

  if (rank == 0) {
    printf("%19s %10s %10s %10s | %10s %10s %10s\n", "",
	   "avg", "min", "max", "cum.time", "cum.cnt", "cum.per");
    printf("%19s %10s %10s %10s | %10s %10s %10s\n", "",
	   "ms", "ms", "ms", "s", "", "ms");
    for (int pr = 0; pr < nr_prof_data; pr++) {
      struct prof_info *pinfo = &prof_globals.info[pr];
      if (pinfo->total_cnt == 0) {
	continue;
      }
      times_avg[pr] /= size;
      
      printf("%-19s %10.2f %10.2f %10.2f | %10.0f %10d %10.2f\n", prof_data[pr].name,
	     times_avg[pr] / 1e3, times_min[pr] / 1e3, times_max[pr] / 1e3,
	     pinfo->total_time / 1e6, pinfo->total_cnt, pinfo->total_time / pinfo->total_cnt / 1e3);
    }
  }

  for (int pr = 0; pr < nr_prof_data; pr++) {
    struct prof_info *pinfo = &prof_globals.info[pr];
    pinfo->time = 0.;
    pinfo->cnt = 0;
  }
}
