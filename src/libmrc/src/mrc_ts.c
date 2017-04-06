
#include <mrc_ts_private.h>
#include <mrc_ts_monitor_private.h>

#include <mrc_ts_monitor.h>
#include <mrc_params.h>
#include <mrc_common.h>
#include <assert.h>

#define mrc_ts_ops(ts) ((struct mrc_ts_ops *) ts->obj.ops)

static void
_mrc_ts_create(struct mrc_ts *ts)
{
  INIT_LIST_HEAD(&ts->monitors);
}

static void
_mrc_ts_setup(struct mrc_ts *ts)
{
  ts->tnorm = ts->norm_time / ts->norm_time_scale;
}

static void
_mrc_ts_view(struct mrc_ts *ts)
{
  if (ts->n == 0)
    return;

  MPI_Comm comm = mrc_ts_comm(ts);
  mpi_printf(comm, "nr_steps      = %d\n", ts->n);
  mpi_printf(comm, "nr_rhsf_evals = %d\n", ts->nr_rhsf_evals);
}

float
mrc_ts_time(struct mrc_ts *ts)
{
  return ts->time;
}

float
mrc_ts_dt(struct mrc_ts *ts)
{
  return ts->dt;
}

int
mrc_ts_step_number(struct mrc_ts *ts)
{
  return ts->n;
}

void
mrc_ts_set_time(struct mrc_ts *ts, float time)
{
  ts->time = time;
}

void
mrc_ts_set_dt(struct mrc_ts *ts, float dt)
{
  ts->dt = dt;
}

void
mrc_ts_set_step_number(struct mrc_ts *ts, int n)
{
  ts->n = n;
}

void
mrc_ts_set_solution(struct mrc_ts *ts, struct mrc_obj *x)
{
  ts->x = x;

  ts->vec_duplicate =
    (struct mrc_obj *(*)(struct mrc_obj *)) mrc_obj_get_method(x, "duplicate");
  assert(ts->vec_duplicate);
  ts->vec_copy =
    (void (*)(struct mrc_obj *, struct mrc_obj *)) mrc_obj_get_method(x, "copy");
  assert(ts->vec_copy);
  ts->vec_axpy =
    (void (*)(struct mrc_obj *, float, struct mrc_obj *)) mrc_obj_get_method(x, "axpy");
  assert(ts->vec_axpy);
  ts->vec_waxpy =
    (void (*)(struct mrc_obj *, float, struct mrc_obj *, struct mrc_obj *)) mrc_obj_get_method(x, "waxpy");
  assert(ts->vec_waxpy);
  ts->vec_norm =
    (float (*)(struct mrc_obj *)) mrc_obj_get_method(x, "norm");
  assert(ts->vec_norm);
  ts->vec_set =
    (void (*)(struct mrc_obj *, float)) mrc_obj_get_method(x, "set");
  assert(ts->vec_set);
}

void
mrc_ts_set_context(struct mrc_ts *ts, struct mrc_obj *ctx_obj)
{
  ts->ctx_obj = ctx_obj;

  mrc_void_func_t rhsf = mrc_obj_get_method(ctx_obj, "rhsf");
  if (rhsf) {
    mrc_ts_set_rhs_function(ts, (void (*)(void *, struct mrc_obj *, float, struct mrc_obj *)) rhsf,
				 ctx_obj);
  }
}

void
mrc_ts_set_rhs_function(struct mrc_ts *ts,
			void (*rhsf)(void *ctx, struct mrc_obj *rhs, float t,
				     struct mrc_obj *x),
			void *ctx)
{
  ts->rhsf = rhsf;
  ts->rhsf_ctx = ctx;
}

void
mrc_ts_set_step_function(struct mrc_ts *ts,
			 void (*stepf)(void *ctx, struct mrc_ts *ts,
				       struct mrc_obj *x),
			 void *ctx)
{
  ts->stepf = stepf;
  ts->stepf_ctx = ctx;
}

void
mrc_ts_set_get_dt_function(struct mrc_ts *ts,
			   double (*get_dt_f)(void *ctx, struct mrc_ts *ts,
					      struct mrc_obj *x),
			   void *ctx)
{
  ts->get_dt_f = get_dt_f;
  ts->get_dt_f_ctx = ctx;
}

void
mrc_ts_set_pre_step_function(struct mrc_ts *ts,
			     void (*pre_step)(void *ctx, struct mrc_ts *ts,
					      struct mrc_obj *x),
			     void *ctx)
{
  ts->pre_step = pre_step;
  ts->pre_step_ctx = ctx;
}

void
mrc_ts_set_post_step_function(struct mrc_ts *ts,
			     void (*post_step)(void *ctx, struct mrc_ts *ts,
					       struct mrc_obj *x),
			     void *ctx)
{
  ts->post_step = post_step;
  ts->post_step_ctx = ctx;
}

void
mrc_ts_add_monitor(struct mrc_ts *ts, struct mrc_ts_monitor *mon)
{
  list_add_tail(&mon->monitors_entry, &ts->monitors);
  mrc_ts_add_child(ts, (struct mrc_obj *) mon);
}

void
mrc_ts_step(struct mrc_ts *ts)
{
  if (ts->pre_step) {
    ts->pre_step(ts->pre_step_ctx, ts, ts->x);
  }

  if (ts->get_dt_f) {
    ts->dt = ts->get_dt_f(ts->rhsf_ctx, ts, ts->x);
    if (ts->time + ts->dt > ts->max_time) {
      ts->dt = ts->max_time - ts->time;
    }
  }

  assert(mrc_ts_ops(ts)->step);
  mrc_ts_ops(ts)->step(ts);

  if (ts->post_step) {
    ts->post_step(ts->post_step_ctx, ts, ts->x);
  }

}

void
mrc_ts_solve(struct mrc_ts *ts)
{
  float _dummyflt = 0.0;
  if (ts->ctx_obj) {
    if (mrc_obj_get_param_float(ts->ctx_obj, "max_time", &_dummyflt) == 0) {
      mrc_obj_set_param_float(ts->ctx_obj, "max_time", ts->max_time);
    }
  }
  
  if (mrc_ts_ops(ts)->solve) {
    mrc_ts_ops(ts)->solve(ts);
    return;
  }

  while ((ts->time < ts->max_time) && (ts->n < ts->max_steps)) {
    mrc_ts_monitors_run(ts);
    mrc_ts_step(ts);
    ts->time += ts->dt;
    ts->n++;
  }

  mrc_ts_monitors_run(ts);
}

// ======================================================================
// helpers for subclasses

void
mrc_ts_monitors_run(struct mrc_ts *ts)
{
  struct mrc_ts_monitor *mon;
  __list_for_each_entry(mon, &ts->monitors, monitors_entry, struct mrc_ts_monitor) {
    mrc_ts_monitor_run(mon, ts);
  }
}

void
mrc_ts_rhsf(struct mrc_ts *ts, struct mrc_obj *rhs, float time,
	    struct mrc_obj *x)
{
  assert(ts->rhsf);
  ts->rhsf(ts->rhsf_ctx, rhs, time, x);
  ts->nr_rhsf_evals++;
}

// ======================================================================
// mrc_ts_create_std
//
// set up a standard mrc_ts with diag and output

struct mrc_ts *
mrc_ts_create_std(MPI_Comm comm,
		  void (*diagf)(void *ctx, float time, struct mrc_obj *x, FILE *file),
		  void *diagf_ctx)
{
  struct mrc_ts *ts = mrc_ts_create(comm);

  struct mrc_ts_monitor *mon_output =
    mrc_ts_monitor_create(mrc_ts_comm(ts));
  mrc_ts_monitor_set_type(mon_output, "output");
  mrc_ts_monitor_set_name(mon_output, "mrc_ts_output");
  mrc_ts_add_monitor(ts, mon_output);

  if (diagf) {
    struct mrc_ts_monitor *mon_diag =
      mrc_ts_monitor_create(mrc_ts_comm(ts));
    mrc_ts_monitor_set_type(mon_diag, "diag");
    mrc_ts_monitor_set_name(mon_diag, "mrc_ts_diag");
    mrc_ts_monitor_diag_set_function(mon_diag, diagf, diagf_ctx);
    mrc_ts_add_monitor(ts, mon_diag);
  }

  return ts;
}

// ======================================================================
// mrc_ts_init

static void
mrc_ts_init()
{
  mrc_class_register_subclass(&mrc_class_mrc_ts, &mrc_ts_step_ops);
  mrc_class_register_subclass(&mrc_class_mrc_ts, &mrc_ts_ode45_ops);
  mrc_class_register_subclass(&mrc_class_mrc_ts, &mrc_ts_rk2_ops);
  mrc_class_register_subclass(&mrc_class_mrc_ts, &mrc_ts_rk4_ops);
  
#ifdef HAVE_PETSC
  mrc_class_register_subclass(&mrc_class_mrc_ts, &mrc_ts_petsc_ops);
#endif

  // Did this break or something? 'Cause, I could fix it, if I had known it was broken. 
  // This is why we have a bug tracking system people...
  //  mrc_class_register_subclass(&mrc_class_mrc_ts, &mrc_ts_rkf45_ops);
}

// ======================================================================
// mrc_ts class

#define VAR(x) (void *)offsetof(struct mrc_ts, x)
static struct param mrc_ts_param_descr[] = {
  { "max_time"       , VAR(max_time)       , PARAM_FLOAT(1.)    },
  { "max_steps"      , VAR(max_steps)      , PARAM_INT(100000)  },
  { "dt"             , VAR(dt)             , PARAM_FLOAT(1e-2)  },
  { "norm_time"      , VAR(norm_time)      , PARAM_DOUBLE(1.)   },
  { "norm_time_scale", VAR(norm_time_scale), PARAM_DOUBLE(1.)   },
  {},
};
#undef VAR

struct mrc_class_mrc_ts mrc_class_mrc_ts = {
  .name         = "mrc_ts",
  .size         = sizeof(struct mrc_ts),
  .param_descr  = mrc_ts_param_descr,
  .init         = mrc_ts_init,
  .create       = _mrc_ts_create,
  .setup        = _mrc_ts_setup,
  .view         = _mrc_ts_view,
};

