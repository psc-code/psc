
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

// ======================================================================

struct mrc_ts_rk4 {
  struct mrc_obj *xt;
  struct mrc_obj *xk[4];
};

static void
mrc_ts_rk4_setup(struct mrc_ts *ts)
{
  struct mrc_ts_rk4 *rk4 = mrc_to_subobj(ts, struct mrc_ts_rk4);
  
  assert(ts->x);
  rk4->xt = mrc_ts_vec_duplicate(ts, ts->x);
  for (int k = 0; k < 4; k++) {
    rk4->xk[k] = mrc_ts_vec_duplicate(ts, ts->x);
  }
}

static void
mrc_ts_rk4_destroy(struct mrc_ts *ts)
{
  struct mrc_ts_rk4 *rk4 = mrc_to_subobj(ts, struct mrc_ts_rk4);

  mrc_obj_destroy(rk4->xt);
  for (int k = 0; k < 4; k++) {
    mrc_obj_destroy(rk4->xk[k]);
  }
}

static void
mrc_ts_rk4_step(struct mrc_ts *ts)
{
  struct mrc_ts_rk4 *rk4 = mrc_to_subobj(ts, struct mrc_ts_rk4);

  struct mrc_obj *x = ts->x;
  struct mrc_obj *xt = rk4->xt;
  struct mrc_obj **xk = rk4->xk;

  mrc_ts_rhsf(ts, xk[0], ts->time, x);

  mrc_ts_vec_waxpy(ts, xt, .5 * ts->dt, xk[0], x);
  mrc_ts_rhsf(ts, xk[1], ts->time + .5 * ts->dt, xt);
    
  mrc_ts_vec_waxpy(ts, xt, .5 * ts->dt, xk[1], x);
  mrc_ts_rhsf(ts, xk[2], ts->time + .5 * ts->dt, xt);
  
  mrc_ts_vec_waxpy(ts, xt, ts->dt, xk[2], x);
  mrc_ts_rhsf(ts, xk[3], ts->time + ts->dt, xt);

  mrc_ts_vec_axpy(ts, x, 1./6. * ts->dt, xk[0]);
  mrc_ts_vec_axpy(ts, x, 1./3. * ts->dt, xk[1]);
  mrc_ts_vec_axpy(ts, x, 1./3. * ts->dt, xk[2]);
  mrc_ts_vec_axpy(ts, x, 1./6. * ts->dt, xk[3]);
}

struct mrc_ts_ops mrc_ts_rk4_ops = {
  .name             = "rk4",
  .size             = sizeof(struct mrc_ts_rk4),
  .setup            = mrc_ts_rk4_setup,
  .destroy          = mrc_ts_rk4_destroy,
  .step             = mrc_ts_rk4_step,
};
