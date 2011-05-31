
#ifndef MRC_TS_MONITOR_PRIVATE_H
#define MRC_TS_MONITOR_PRIVATE_H

#include <mrc_ts_monitor.h>

// ======================================================================
// mrc_ts_monitor

struct mrc_ts_monitor {
  struct mrc_obj obj;
  // parameters
  int every_steps;
  float every_time;

  int next_step;
  float next_time;
  list_t monitors_entry;
};

struct mrc_ts_monitor_ops {
  MRC_SUBCLASS_OPS(struct mrc_ts_monitor);
  void (*run)(struct mrc_ts_monitor *, struct mrc_ts *);
};

extern struct mrc_ts_monitor_ops mrc_ts_monitor_output_ops;
extern struct mrc_ts_monitor_ops mrc_ts_monitor_diag_ops;

#endif
