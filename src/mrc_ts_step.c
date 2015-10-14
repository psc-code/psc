
#include <mrc_ts_private.h>

#include <assert.h>

// ======================================================================
// mrc_ts subclass "step"
//
// the simplest possible timestepper, because it doesn't actually do
// anything other than delegate to a function the user has set.

// ----------------------------------------------------------------------
// mrc_ts_step_step

static void
mrc_ts_step_step(struct mrc_ts *ts)
{
  assert(ts->stepf);
  ts->stepf(ts->rhsf_ctx, ts, ts->x);

  // FIXME, this should go -> mrc_ts_step(), and be optional
  mpi_printf(mrc_ts_comm(ts), "step=%i time=%g dt=%e\n",
	     ts->n + 1, ts->time + ts->dt, ts->dt);
}

// ----------------------------------------------------------------------
// mrc_ts_step_ops

struct mrc_ts_ops mrc_ts_step_ops = {
  .name             = "step",
  .step             = mrc_ts_step_step,
};

