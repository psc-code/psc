
#include <mrc_ts_private.h>
#include <mrc_ts_monitor_private.h>

#include <mrc_ts_monitor.h>
#include <mrc_params.h>
#include <assert.h>

#define mrc_ts_ops(ts) ((struct mrc_ts_ops *) ts->obj.ops)

static void
_mrc_ts_create(struct mrc_ts *ts)
{
  INIT_LIST_HEAD(&ts->monitors);
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

void
mrc_ts_set_dt(struct mrc_ts *ts, float dt)
{
  ts->dt = dt;
}

void
mrc_ts_set_solution(struct mrc_ts *ts, struct mrc_f1 *x)
{
  ts->x = x;
}

void
mrc_ts_set_rhs_function(struct mrc_ts *ts,
			void (*rhsf)(void *ctx, struct mrc_f1 *x,
				     struct mrc_f1 *rhs),
			void *ctx)
{
  ts->rhsf = rhsf;
  ts->ctx = ctx;
}

void
mrc_ts_add_monitor(struct mrc_ts *ts, struct mrc_ts_monitor *mon)
{
  list_add_tail(&mon->monitors_entry, &ts->monitors);
  mrc_ts_add_child(ts, (struct mrc_obj *) mon);
}

static void
mrc_ts_step(struct mrc_ts *ts)
{
  assert(mrc_ts_ops(ts)->step);
  mrc_ts_ops(ts)->step(ts);
}

void
mrc_ts_solve(struct mrc_ts *ts)
{
  if (mrc_ts_ops(ts)->solve) {
    mrc_ts_ops(ts)->solve(ts);
    return;
  }

  while ((ts->time < ts->max_time) && (ts->n < ts->max_steps)) {
    if (ts->time + ts->dt > ts->max_time) {
      ts->dt = ts->max_time - ts->time;
    }

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
  list_for_each_entry(mon, &ts->monitors, monitors_entry) {
    mrc_ts_monitor_run(mon, ts);
  }
}

void
mrc_ts_rhsf(struct mrc_ts *ts, struct mrc_f1 *rhs, float time,
	    struct mrc_f1 *x)
{
  ts->rhsf(ts->ctx, rhs, x);
  ts->nr_rhsf_evals++;
}

// ======================================================================
// mrc_ts_create_std
//
// set up a standard mrc_ts with diag and output

struct mrc_ts *
mrc_ts_create_std(MPI_Comm comm,
		  void (*diagf)(void *ctx, float time, struct mrc_f1 *x, FILE *file))
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
    mrc_ts_monitor_diag_set_function(mon_diag, diagf);
    mrc_ts_add_monitor(ts, mon_diag);
  }

  return ts;
}

// ======================================================================
// mrc_ts_init

static void
mrc_ts_init()
{
  mrc_class_register_subclass(&mrc_class_mrc_ts, &mrc_ts_ode45_ops);
  mrc_class_register_subclass(&mrc_class_mrc_ts, &mrc_ts_rk2_ops);
  mrc_class_register_subclass(&mrc_class_mrc_ts, &mrc_ts_rk4_ops);
}

// ======================================================================
// mrc_ts class

#define VAR(x) (void *)offsetof(struct mrc_ts, x)
static struct param mrc_ts_param_descr[] = {
  { "max_time"      , VAR(max_time)      , PARAM_FLOAT(1.)           },
  { "max_steps"     , VAR(max_steps)     , PARAM_INT(100000)         },
  { "dt"            , VAR(dt)            , PARAM_FLOAT(1e-2)           },
  {},
};
#undef VAR

struct mrc_class_mrc_ts mrc_class_mrc_ts = {
  .name         = "mrc_ts",
  .size         = sizeof(struct mrc_ts),
  .param_descr  = mrc_ts_param_descr,
  .init         = mrc_ts_init,
  .create       = _mrc_ts_create,
  .view         = _mrc_ts_view,
};

