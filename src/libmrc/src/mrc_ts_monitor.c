
#include <mrc_ts_monitor_private.h>

#include <mrc_ts_private.h> // FIXME?
#include <mrc_params.h>
#include <assert.h>

// ======================================================================

#define mrc_ts_monitor_ops(mon) ((struct mrc_ts_monitor_ops *) mon->obj.ops)

void
mrc_ts_monitor_run(struct mrc_ts_monitor *mon, struct mrc_ts *ts)
{
  assert(mrc_ts_monitor_ops(mon)->run);

  double time = ts->time * ts->tnorm;
  if ((mon->every_steps > 0 && ts->n >= mon->next_step) ||
      (mon->every_time > 0. && time >= mon->next_time)) {
    mrc_ts_monitor_ops(mon)->run(mon, ts);
    mon->next_step = ts->n + mon->every_steps;
    while (mon->next_time <= time) {
      mon->next_time += mon->every_time;
    }
  } else if (ts->time >= ts->max_time) {
    // do final output at the end
    mrc_ts_monitor_ops(mon)->run(mon, ts);
  }
}

// ======================================================================
// mrc_ts_monitor_init

static void
mrc_ts_monitor_init()
{
  mrc_class_register_subclass(&mrc_class_mrc_ts_monitor, &mrc_ts_monitor_output_ops);
  mrc_class_register_subclass(&mrc_class_mrc_ts_monitor, &mrc_ts_monitor_diag_ops);
}

// ======================================================================
// mrc_ts_monitor class

#define VAR(x) (void *)offsetof(struct mrc_ts_monitor, x)
static struct param mrc_ts_monitor_descr[] = {
  { "every_steps"    , VAR(every_steps)    , PARAM_INT(0)             },
  { "every_time"     , VAR(every_time)     , PARAM_FLOAT(0.1)         },
  {},
};
#undef VAR

struct mrc_class_mrc_ts_monitor mrc_class_mrc_ts_monitor = {
  .name         = "mrc_ts_monitor",
  .size         = sizeof(struct mrc_ts_monitor),
  .param_descr  = mrc_ts_monitor_descr,
  .init         = mrc_ts_monitor_init,
};

