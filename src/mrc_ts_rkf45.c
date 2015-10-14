

#include <mrc_ts_private.h>
#include <math.h>

#include <assert.h>

// Runge-Kutta Fehlburg 4th/5th order adaptive integrator
struct mrc_ts_rkf45 {
  double max_error; // absolute error for now
  double dt_min; // minimum stepsize before killing
  struct mrc_obj *xt;
  struct mrc_obj *xt2; // need 2 temp variables for adaptive
  struct mrc_obj *xk[6];
};

#define VAR(x) (void *)offsetof(struct mrc_ts_rkf45, x)
static struct param mrc_ts_rkf45_param_descr[] = {
  { "max_error"      , VAR(max_error)      , PARAM_DOUBLE(1e-7)           },
  { "dt_min"         , VAR(dt_min)         , PARAM_DOUBLE(1e-10)           },
  {},
};
#undef VAR

static void
mrc_ts_rkf45_setup(struct mrc_ts *ts)
{
  struct mrc_ts_rkf45 *rkf45 = mrc_to_subobj(ts, struct mrc_ts_rkf45);
  
  assert(ts->x);
  rkf45->xt = mrc_ts_vec_duplicate(ts, ts->x);
  rkf45->xt2 = mrc_ts_vec_duplicate(ts, ts->x);
  for (int k = 0; k < 6; k++) {
    rkf45->xk[k] = mrc_ts_vec_duplicate(ts, ts->x);
  }

  mrc_ts_setup_super(ts);
}

static void
mrc_ts_rkf45_destroy(struct mrc_ts *ts)
{
  struct mrc_ts_rkf45 *rkf45 = mrc_to_subobj(ts, struct mrc_ts_rkf45);

  mrc_obj_destroy(rkf45->xt);
  mrc_obj_destroy(rkf45->xt2);
  for (int k = 0; k < 6; k++) {
    mrc_obj_destroy(rkf45->xk[k]);
  }
}

static void
mrc_ts_rkf45_step(struct mrc_ts *ts)
{
  struct mrc_ts_rkf45 *rkf45 = mrc_to_subobj(ts, struct mrc_ts_rkf45);

  struct mrc_obj *x = ts->x;
  struct mrc_obj *xt = rkf45->xt;
  struct mrc_obj *x_4th = rkf45->xt2;
  struct mrc_obj **xk = rkf45->xk;

  // Bring in coefficients
#include "rkf45_butcher.h"
  // Create arrays of  scalars for the coefficients
  double first[1]  = { A10 };
  double second[2] = { A20 , A21 };
  double third[3]  = { A30 , A31 , A32 };
  double fourth[4] = { A40 , A41 , A42 , A43 };
  double fifth[5]  = { A50 , A51, A52 , A53 , A54 };

  double *coeff[5] = { first, second, third, fourth, fifth };
  double time_coeff[5] = { C1 , C2 , C3 , C4 , C5 };
  // Fourth order coeff
  double coeff4th[6] = { B4_0 , B4_1 , B4_2 , B4_3 , B4_4 , B4_5 };
  // Fifth order coeff
  double coeff5th[6] = { B5_0 , B5_1 , B5_2 , B5_3 , B5_4 , B5_5 };

  // zero order
  mrc_ts_rhsf(ts, xk[0], ts->time, x);
  
  // Now we can loop to fill out the k values
  // Also, these index names a little goofy.
  for (int ik = 0; ik < 5; ik++) {
    // using a waxpy here means we don't need to zero xt
    mrc_ts_vec_waxpy(ts, xt, ts->dt * coeff[ik][0], xk[0], x);
    // dd all the lower k's times their coeffficients
    for (int n = 1; n <= ik; n++) {
      mrc_ts_vec_axpy(ts, xt, ts->dt * coeff[ik][n], xk[n]);
    }
    mrc_ts_rhsf(ts, xk[ik+1], ts->time + time_coeff[ik] * ts->dt, xt);
  }

  // Calculate the sum for 4th order solution,
  // which will be stored in x_4th,
  // and which we will add to x for the final solution
  mrc_ts_vec_set(ts, x_4th, 0.0);
  
  for (int i = 0; i < 6; i++) {
    if (coeff4th[i] != 0.0) {
      mrc_ts_vec_axpy(ts, x_4th, coeff4th[i], xk[i]);
    }
  }
  
  // Rather than actually calculating the 5th order solution, 
  // which we only use to check the error,
  // we'll take the difference of the 4th and fifth order
  // and store into rv
  
  // Just sum and store the first k, so we don't have to
  // zero rv
  mrc_ts_vec_waxpy(ts, xt, -1. * coeff5th[0], xk[0], x_4th);
  for (int i = 1; i < 6; i++) {
    if (coeff5th[i] != 0.0) {
      mrc_ts_vec_axpy(ts, xt, -1.* coeff5th[i], xk[i]);
    }
  }
  
  // calc error and adjust dt
  float norm;
  norm = ts->dt * mrc_ts_vec_norm(ts, xt);
  double dt_fac = exp(log((rkf45->max_error) / norm) / ((6) + 1)) * 0.9;

  if (dt_fac > 5.0) dt_fac = 5.0;

  double old_dt = ts->dt;
  ts->dt *= dt_fac;
  
  assert( ts->dt > rkf45->dt_min );

  // mprintf("t %d:%g dt = %g norm=%g\n", ts->n, ts->time, ts->dt, norm);
  float tol = rkf45->max_error * fmaxf(mrc_ts_vec_norm(ts, x), 1.); // allowable error 

  if ( norm <= tol) 
    { 
      mrc_ts_vec_axpy(ts, x, old_dt, x_4th);
      // need to hack, because time gets incremented
      // outside of single_step, but we want to set the new dt here
      ts->time += old_dt - ts->dt;
    } else {
    // if we don't converge, call recursively until we do a good job,
    // or hit dt_min
    mrc_ts_rkf45_step(ts);
  }
}

struct mrc_ts_ops mrc_ts_rkf45_ops = {
  .name             = "rkf45",
  .size             = sizeof(struct mrc_ts_rkf45),
  .param_descr      = mrc_ts_rkf45_param_descr,
  .setup            = mrc_ts_rkf45_setup,
  .destroy          = mrc_ts_rkf45_destroy,
  .step             = mrc_ts_rkf45_step,
};
