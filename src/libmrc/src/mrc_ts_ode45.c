
#include <mrc_ts_private.h>

#include <mrc_params.h>
#include <math.h>
#include <assert.h>

// Adapted from octave's ode45.m

struct mrc_ts_ode45 {
  float tol;
  float safety_factor;
  float dt_min;
  float dt_max;
  float pow;

  int nr_rejected_steps;

  struct mrc_obj *xk[7];
  struct mrc_obj *x4;
  struct mrc_obj *x5;
  struct mrc_obj *gamma1;
};

#define VAR(x) (void *)offsetof(struct mrc_ts_ode45, x)
static struct param mrc_ts_ode45_descr[] = {
  { "tol"             , VAR(tol)           , PARAM_FLOAT(1e-6)      },
  { "safety_factor"   , VAR(safety_factor) , PARAM_FLOAT(.8)        },
  { "dt_min"          , VAR(dt_min)        , PARAM_FLOAT(0.)        },
  { "dt_max"          , VAR(dt_max)        , PARAM_FLOAT(0.)        },
  {},
};
#undef VAR

// ======================================================================

// Dormand-Prince 4(5) coefficients:
static const float a[7][6] = {
  { },
  { 1./5, },
  { 3./40., 9./40. },
  { 44./45., -56./15., 32./9. },
  { 19372./6561., -25360./2187., 64448./6561., -212./729. },
  { 9017./3168., -355./33., 46732./5247., 49./176., -5103./18656. },
  { 35./384., 0., 500./1113., 125./192., -2187./6784., 11./84. }
};

// 4th order b-coefficients
static const float b4[7] = {
  5179./57600., 0., 7571./16695., 393./640., -92097./339200., 187./2100., 1./40.
};

// 5th order b-coefficients
static const float b5[7] = {
  35./384., 0., 500./1113., 125./192., -2187./6784., 11./84., 0.
};

static float c[7];

// ----------------------------------------------------------------------

static void
mrc_ts_ode45_setup(struct mrc_ts *ts)
{
  // calculating these once would be enough
  for (int i = 0; i < 7; i++) {
    c[i] = 0.;
    for (int j = 0; j < 6; j++) {
      c[i] += a[i][j];
    }
  }

  struct mrc_ts_ode45 *ode45 = mrc_to_subobj(ts, struct mrc_ts_ode45);
  assert(ts->x);
  
  for (int i = 0; i < 7; i++) {
    ode45->xk[i] = mrc_ts_vec_duplicate(ts, ts->x);
  }
  ode45->x4 = mrc_ts_vec_duplicate(ts, ts->x);
  ode45->x5 = mrc_ts_vec_duplicate(ts, ts->x);
  ode45->gamma1 = mrc_ts_vec_duplicate(ts, ts->x);

  ode45->pow = 1./6.; // see p.91 in the Ascher & Petzold reference for more infomation.

  if (ode45->dt_max == 0.) {
    ode45->dt_max = (ts->max_time - ts->time) / 2.5;
  }
  if (ode45->dt_min == 0.) {
    ode45->dt_min = (ts->max_time - ts->time) / 1e9;
  }
  if (ts->dt == 0.) {
    ts->dt = (ts->max_time - ts->time) / 100.; // initial guess at a step size
  }

  mrc_ts_setup_super(ts);
}

static void
mrc_ts_ode45_destroy(struct mrc_ts *ts)
{
  struct mrc_ts_ode45 *ode45 = mrc_to_subobj(ts, struct mrc_ts_ode45);

  for (int i = 0; i < 7; i++) {
    mrc_obj_destroy(ode45->xk[i]);
  }
  mrc_obj_destroy(ode45->x4);
  mrc_obj_destroy(ode45->x5);
  mrc_obj_destroy(ode45->gamma1);
}

static void
mrc_ts_ode45_view(struct mrc_ts *ts)
{
  struct mrc_ts_ode45 *ode45 = mrc_to_subobj(ts, struct mrc_ts_ode45);

  if (ts->n == 0)
    return;

  MPI_Comm comm = mrc_ts_comm(ts);
  mpi_printf(comm, "nr_rejected_steps = %d\n", ode45->nr_rejected_steps);
}

static void
mrc_ts_ode45_solve(struct mrc_ts *ts)
{
  struct mrc_ts_ode45 *ode45 = mrc_to_subobj(ts, struct mrc_ts_ode45);
  struct mrc_obj *x = ts->x;
  struct mrc_obj **xk = ode45->xk;
  struct mrc_obj *x4 = ode45->x4;
  struct mrc_obj *x5 = ode45->x5;
  struct mrc_obj *gamma1 = ode45->gamma1;

  mrc_ts_rhsf(ts, xk[0], ts->time, x);

  while ((ts->time < ts->max_time) && (ts->dt >= ode45->dt_min)) {
    if (ts->time + ts->dt > ts->max_time) {
      ts->dt = ts->max_time - ts->time;
    }

    mrc_ts_monitors_run(ts);

    // Compute the slopes by computing the k(:,j+1)'th column based on 
    // the previous k(:,1:j) columns
    // notes: k_ needs to end up as an Nxs, a_ is 7x6, which is s by (s-1),
    // s is the number of intermediate RK stages on [t (t+h)] (Dormand-Prince has s=7 stages)
    for (int j = 0; j < 6; j++) {
      // k_(:,j+1) = feval(FUN, t+c_(j+1)*h, x+h*k_(:,1:j)*a_(j+1,1:j) );
      mrc_ts_vec_copy(ts, gamma1, x);
      for (int k = 0; k <= j; k++) {
	mrc_ts_vec_axpy(ts, gamma1, ts->dt * a[j+1][k], xk[k]);
      }
      float time = ts->time + c[j+1] * ts->dt;
      mrc_ts_rhsf(ts, xk[j+1], time, gamma1);
    }

    // compute the 4th order estimate
    //x4=x + h* (k_*b4_);
    mrc_ts_vec_copy(ts, x4, x);
    for (int k = 0; k < 7; k++) {
      mrc_ts_vec_axpy(ts, x4, ts->dt * b4[k], xk[k]);
    }

    // compute the 5th order estimate
    //x5=x + h*(k_*b5_);
    mrc_ts_vec_copy(ts, x5, x);
    for (int k = 0; k < 7; k++) {
      mrc_ts_vec_axpy(ts, x5, ts->dt * b5[k], xk[k]);
    }

    // estimate the local truncation error
    mrc_ts_vec_waxpy(ts, gamma1, -1., x4, x5);

    // Estimate the error and the acceptable error
    float delta = mrc_ts_vec_norm(ts, gamma1); // actual error
    float tau = ode45->tol * fmaxf(mrc_ts_vec_norm(ts, x), 1.); // allowable error

    // Update the solution only if the error is acceptable
    if (delta <= tau) {
      ts->time += ts->dt;
      ts->n++;
      // using the higher order estimate is called 'local extrapolation'
      mrc_ts_vec_copy(ts, x, x5);
    } else {
      ode45->nr_rejected_steps++;
    }

    // Update the step size
    if (delta == 0.) {
      delta = 1e-16;
    }
    ts->dt = fminf(ode45->dt_max, 
		   ode45->safety_factor * ts->dt * pow(tau/delta, ode45->pow));

    // Assign the last stage for x(k) as the first stage for computing x(k+1).
    // This is part of the Dormand-Prince pair caveat.
    // k_(:,7) has already been computed, so use it instead of recomputing it
    // again as k_(:,1) during the next step.
    mrc_ts_vec_copy(ts, xk[0], xk[6]);

    //    mprintf("t %d:%g dt = %g\n", ts->n, ts->time, ts->dt);
  }

  mrc_ts_monitors_run(ts);
}

// ----------------------------------------------------------------------=
// mrc_ts_ode45

struct mrc_ts_ops mrc_ts_ode45_ops = {
  .name             = "ode45",
  .size             = sizeof(struct mrc_ts_ode45),
  .param_descr      = mrc_ts_ode45_descr,
  .setup            = mrc_ts_ode45_setup,
  .destroy          = mrc_ts_ode45_destroy,
  .view             = mrc_ts_ode45_view,
  .solve            = mrc_ts_ode45_solve,
};
