
#include <mrc_ts_monitor_private.h>

#include <mrc_ts_private.h> // FIXME?
#include <mrc_io.h>

struct mrc_ts_monitor_output {
  struct mrc_io *io;
  int nr;
};

static void
mrc_ts_monitor_output_create(struct mrc_ts_monitor *mon)
{
  struct mrc_ts_monitor_output *out =
    mrc_to_subobj(mon, struct mrc_ts_monitor_output);

  out->io = mrc_io_create(mrc_ts_monitor_comm(mon));
  mrc_ts_monitor_add_child(mon, (struct mrc_obj *) out->io);
}

static void
mrc_ts_monitor_output_run(struct mrc_ts_monitor *mon, struct mrc_ts *ts)
{
  struct mrc_ts_monitor_output *out =
    mrc_to_subobj(mon, struct mrc_ts_monitor_output);

  mpi_printf(mrc_ts_monitor_comm(mon), "Writing output %d (time = %g)\n",
	     out->nr, ts->time);
  mrc_io_open(out->io, "w", out->nr, ts->time);
  mrc_obj_write(ts->x, out->io);
  mrc_io_close(out->io);
  out->nr++;
}

struct mrc_ts_monitor_ops mrc_ts_monitor_output_ops = {
  .name             = "output",
  .size             = sizeof(struct mrc_ts_monitor_output),
  .create           = mrc_ts_monitor_output_create,
  .run              = mrc_ts_monitor_output_run,
};

