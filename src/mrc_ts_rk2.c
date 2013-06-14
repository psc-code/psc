
#include <mrc_ts_private.h>

#include <assert.h>

struct mrc_ts_rk2 {
  struct mrc_obj *rhs;
  struct mrc_obj *xm;
};

static void
mrc_ts_rk2_setup(struct mrc_ts *ts)
{
  struct mrc_ts_rk2 *rk2 = mrc_to_subobj(ts, struct mrc_ts_rk2);
  
  assert(ts->x);
  rk2->xm = mrc_ts_vec_duplicate(ts, ts->x);
  rk2->rhs = mrc_ts_vec_duplicate(ts, ts->x);

  mrc_ts_setup_super(ts);
}

static void
mrc_ts_rk2_destroy(struct mrc_ts *ts)
{
  struct mrc_ts_rk2 *rk2 = mrc_to_subobj(ts, struct mrc_ts_rk2);

  mrc_obj_destroy(rk2->xm);
  mrc_obj_destroy(rk2->rhs);
}

static void
mrc_ts_rk2_step(struct mrc_ts *ts)
{
  struct mrc_ts_rk2 *rk2 = mrc_to_subobj(ts, struct mrc_ts_rk2);

  struct mrc_obj *x = ts->x;
  struct mrc_obj *rhs = rk2->rhs;
  struct mrc_obj *xm = rk2->xm;

  mrc_ts_rhsf(ts, rhs, ts->time, x);
  mrc_ts_vec_waxpy(ts, xm, .5 * ts->dt, rhs, x);
  
  mrc_ts_rhsf(ts, rhs, ts->time + .5 * ts->dt, xm);
  mrc_ts_vec_axpy(ts, x, ts->dt, rhs);
}

struct mrc_ts_ops mrc_ts_rk2_ops = {
  .name             = "rk2",
  .size             = sizeof(struct mrc_ts_rk2),
  .setup            = mrc_ts_rk2_setup,
  .destroy          = mrc_ts_rk2_destroy,
  .step             = mrc_ts_rk2_step,
};

