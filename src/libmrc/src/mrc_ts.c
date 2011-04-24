
#include <mrc_ts_private.h>

#include <mrc_params.h>

static void
_mrc_ts_create(struct mrc_ts *ts)
{
  ts->io = mrc_io_create(mrc_ts_comm(ts));
  mrc_ts_add_child(ts, (struct mrc_obj *) ts->io);
}

static void
_mrc_ts_setup(struct mrc_ts *ts)
{
  ts->f_diag = fopen(ts->diag_filename, "w");

  // RK2
  ts->rhs = mrc_f1_duplicate(ts->x);
  ts->xm = mrc_f1_duplicate(ts->x);
}

static void
_mrc_ts_destroy(struct mrc_ts *ts)
{
  if (ts->f_diag) {
    fclose(ts->f_diag);
  }

  // RK2
  mrc_f1_destroy(ts->rhs);
  mrc_f1_destroy(ts->xm);
}

void
mrc_ts_set_dt(struct mrc_ts *ts, float dt)
{
  ts->dt = dt;
}

void
mrc_ts_set_state(struct mrc_ts *ts, struct mrc_f1 *x)
{
  ts->x = x;
}

void
mrc_ts_set_context(struct mrc_ts *ts, void *ctx)
{
  ts->ctx = ctx;
}

void
mrc_ts_set_rhs_function(struct mrc_ts *ts,
			void (*rhsf)(void *ctx, struct mrc_f1 *x,
				     struct mrc_f1 *rhs))
{
  ts->rhsf = rhsf;
}

void
mrc_ts_set_diag_function(struct mrc_ts *ts,
			 void (*diagf)(void *ctx, float time, struct mrc_f1 *x,
				       FILE *file))
{
  ts->diagf = diagf;
}

static void
mrc_ts_diag(struct mrc_ts *ts, float time)
{
  if (ts->n % ts->diag_every == 0) {
    ts->diagf(ts->ctx, time, ts->x, ts->f_diag);
  }
}

static void
mrc_ts_output(struct mrc_ts *ts, float time)
{
  if (ts->n % ts->out_every == 0) {
    mrc_io_open(ts->io, "w", ts->n / ts->out_every, time);
    mrc_f1_write(ts->x, ts->io);
    mrc_io_close(ts->io);
  }
}

static void
mrc_ts_step(struct mrc_ts *ts)
{
  struct mrc_f1 *x = ts->x;
  struct mrc_f1 *xm = ts->xm;
  struct mrc_f1 *rhs = ts->rhs;

  ts->rhsf(ts->ctx, rhs, x);
  mrc_f1_waxpy(xm, .5 * ts->dt, rhs, x);
  
  ts->rhsf(ts->ctx, rhs, xm);
  mrc_f1_axpy(x, ts->dt, rhs);
}

void
mrc_ts_solve(struct mrc_ts *ts)
{
  for (ts->n = 0; ts->n <= ts->max_steps; ts->n++) {
    mrc_ts_output(ts, ts->n * ts->dt);
    mrc_ts_diag(ts, ts->n * ts->dt);
    mrc_ts_step(ts);
  }
}

#define VAR(x) (void *)offsetof(struct mrc_ts, x)
static struct param mrc_ts_param_descr[] = {
  { "max_steps"     , VAR(max_steps)     , PARAM_INT(10000)          },
  { "out_every"     , VAR(out_every)     , PARAM_INT(1000)           },
  { "diag_every"    , VAR(diag_every)    , PARAM_INT(10)             },
  { "diag_filename" , VAR(diag_filename) , PARAM_STRING("diag.asc")  },
  {},
};
#undef VAR

struct mrc_class_mrc_ts mrc_class_mrc_ts = {
  .name         = "mrc_ts",
  .size         = sizeof(struct mrc_ts),
  .param_descr  = mrc_ts_param_descr,
  .create       = _mrc_ts_create,
  .setup        = _mrc_ts_setup,
  .destroy      = _mrc_ts_destroy,
};

