
#ifndef MRC_TS_PRIVATE_H
#define MRC_TS_PRIVATE_H

#include <mrc_ts.h>

#include <stdio.h>

struct mrc_ts {
  struct mrc_obj obj;
  // parameters
  int max_steps;
  float max_time;

  int n; // current timestep number
  float time; // current integration time
  float dt; // current dt
  struct mrc_f1 *x; // current state vector
  void *ctx; // FIXME, should be mrc_obj?
  void (*rhsf)(void *ctx, struct mrc_f1 *x, struct mrc_f1 *rhs);

  list_t monitors;
  void (*diagf)(void *ctx, float time, struct mrc_f1 *x, FILE *file);

  struct mrc_f1 *rhs;
  int nr_rhsf_evals; // statistics
};

// helpers for subclasses

void mrc_ts_rhsf(struct mrc_ts *ts, struct mrc_f1 *rhs, float time,
		 struct mrc_f1 *x);
void mrc_ts_monitors_run(struct mrc_ts *ts);

struct mrc_ts_ops {
  MRC_SUBCLASS_OPS(struct mrc_ts);
  void (*step)(struct mrc_ts *);
  void (*solve)(struct mrc_ts *);
};

extern struct mrc_ts_ops mrc_ts_rk2_ops;
extern struct mrc_ts_ops mrc_ts_ode45_ops;

#endif
