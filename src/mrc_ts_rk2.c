
#include <mrc_ts_private.h>

struct mrc_ts_rk2 {
  struct mrc_f1 *xm;
};

static void
mrc_ts_rk2_setup(struct mrc_ts *ts)
{
  struct mrc_ts_rk2 *rk2 = mrc_to_subobj(ts, struct mrc_ts_rk2);
  
  rk2->xm = mrc_f1_duplicate(ts->x);
}

static void
mrc_ts_rk2_destroy(struct mrc_ts *ts)
{
  struct mrc_ts_rk2 *rk2 = mrc_to_subobj(ts, struct mrc_ts_rk2);

  mrc_f1_destroy(rk2->xm);
}

static void
mrc_ts_rk2_step(struct mrc_ts *ts)
{
  struct mrc_ts_rk2 *rk2 = mrc_to_subobj(ts, struct mrc_ts_rk2);

  struct mrc_f1 *x = ts->x;
  struct mrc_f1 *rhs = ts->rhs;
  struct mrc_f1 *xm = rk2->xm;

  mrc_ts_rhsf(ts, rhs, ts->time, x);
  mrc_f1_waxpy(xm, .5 * ts->dt, rhs, x);
  
  mrc_ts_rhsf(ts, rhs, ts->time + .5 * ts->dt, xm);
  mrc_f1_axpy(x, ts->dt, rhs);
}

struct mrc_ts_ops mrc_ts_rk2_ops = {
  .name             = "rk2",
  .size             = sizeof(struct mrc_ts_rk2),
  .setup            = mrc_ts_rk2_setup,
  .destroy          = mrc_ts_rk2_destroy,
  .step             = mrc_ts_rk2_step,
};
