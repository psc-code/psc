
#include <mrc_ts_monitor_private.h>

#include <mrc_ts_private.h> // FIXME?
#include <mrc_params.h>

struct mrc_ts_monitor_diag {
  char *filename;

  FILE *file;
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
  struct mrc_ts_monitor_diag *diag =
    mrc_to_subobj(mon, struct mrc_ts_monitor_diag);

  diag->file = fopen(diag->filename, "w");
}

static void
mrc_ts_monitor_diag_destroy(struct mrc_ts_monitor *mon)
{
  struct mrc_ts_monitor_diag *diag =
    mrc_to_subobj(mon, struct mrc_ts_monitor_diag);

  if (diag->file) {
    fclose(diag->file);
  }
}

static void
mrc_ts_monitor_diag_run(struct mrc_ts_monitor *mon, struct mrc_ts *ts)
{
  struct mrc_ts_monitor_diag *diag =
    mrc_to_subobj(mon, struct mrc_ts_monitor_diag);

  mpi_printf(mrc_ts_monitor_comm(mon), "writing diag %d (%g)\n", ts->n, ts->time);
  ts->diagf(ts->ctx, ts->time, ts->x, diag->file);
}

struct mrc_ts_monitor_ops mrc_ts_monitor_diag_ops = {
  .name             = "diag",
  .size             = sizeof(struct mrc_ts_monitor_diag),
  .param_descr      = mrc_ts_monitor_diag_descr,
  .setup            = mrc_ts_monitor_diag_setup,
  .destroy          = mrc_ts_monitor_diag_destroy,
  .run              = mrc_ts_monitor_diag_run,
};

