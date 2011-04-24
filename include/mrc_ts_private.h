
#ifndef MRC_TS_PRIVATE_H
#define MRC_TS_PRIVATE_H

#include <mrc_ts.h>

#include <mrc_io.h>
#include <stdio.h>

struct mrc_ts {
  struct mrc_obj obj;
  // parameters
  int out_every;
  int diag_every;
  int max_steps;
  float max_time;
  char *diag_filename;

  struct mrc_io *io; // for writing state output
  int n; // current timestep number
  float time; // current integration time
  float dt; // current dt
  struct mrc_f1 *x; // current state vector
  void *ctx; // FIXME, should be mrc_obj?
  void (*rhsf)(void *ctx, struct mrc_f1 *x, struct mrc_f1 *rhs);
  FILE *f_diag;
  void (*diagf)(void *ctx, float time, struct mrc_f1 *x, FILE *file);

  struct mrc_f1 *rhs;
};

void mrc_ts_diag(struct mrc_ts *ts);
void mrc_ts_output(struct mrc_ts *ts);

struct mrc_ts_ops {
  MRC_SUBCLASS_OPS(struct mrc_ts);
  void (*step)(struct mrc_ts *);
  void (*solve)(struct mrc_ts *);
};

extern struct mrc_ts_ops mrc_ts_rk2_ops;
extern struct mrc_ts_ops mrc_ts_ode45_ops;

#endif
