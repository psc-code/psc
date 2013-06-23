
#include <mrc_ts_monitor_private.h>

#include <mrc_ts_private.h> // FIXME?
#include <mrc_params.h>

#define to_diag(mon) mrc_to_subobj(mon, struct mrc_ts_monitor_diag)

struct mrc_ts_monitor_diag {
  char *filename;

  FILE *file;
  void (*diagf)(void *ctx, float time, struct mrc_obj *x, FILE *file);
  void *diagf_ctx;
};

#define VAR(x) (void *)offsetof(struct mrc_ts_monitor_diag, x)
static struct param mrc_ts_monitor_diag_descr[] = {
  { "filename"      , VAR(filename)      , PARAM_STRING("diag.asc")             },
  {},
};
#undef VAR

static void
mrc_ts_monitor_diag_setup(struct mrc_ts_monitor *mon)
{
  struct mrc_ts_monitor_diag *diag = to_diag(mon);

  int rank;
  MPI_Comm_rank(mrc_ts_monitor_comm(mon), &rank);
  if (rank == 0) {
    diag->file = fopen(diag->filename, "w");
  }
}

static void
mrc_ts_monitor_diag_destroy(struct mrc_ts_monitor *mon)
{
  struct mrc_ts_monitor_diag *diag = to_diag(mon);

  if (diag->file) {
    fclose(diag->file);
  }
}

void
mrc_ts_monitor_diag_set_function(struct mrc_ts_monitor *mon,
				 void (*diagf)(void *ctx, float time, struct mrc_obj *x,
					       FILE *file),
				 void *diagf_ctx)
{
  struct mrc_ts_monitor_diag *diag = to_diag(mon);

  diag->diagf = diagf;
  diag->diagf_ctx = diagf_ctx;
}

static void
mrc_ts_monitor_diag_run(struct mrc_ts_monitor *mon, struct mrc_ts *ts)
{
  struct mrc_ts_monitor_diag *diag = to_diag(mon);

  mpi_printf(mrc_ts_monitor_comm(mon), "writing diag %d (%g)\n", ts->n, ts->time);
  diag->diagf(diag->diagf_ctx, ts->time, ts->x, diag->file);
  fflush(diag->file);
}

struct mrc_ts_monitor_ops mrc_ts_monitor_diag_ops = {
  .name             = "diag",
  .size             = sizeof(struct mrc_ts_monitor_diag),
  .param_descr      = mrc_ts_monitor_diag_descr,
  .setup            = mrc_ts_monitor_diag_setup,
  .destroy          = mrc_ts_monitor_diag_destroy,
  .run              = mrc_ts_monitor_diag_run,
};

