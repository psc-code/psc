
#ifndef MRC_TS_H
#define MRC_TS_H

#include <mrc_obj.h>

#include <mrc_fld.h>
#include <stdio.h>

struct mrc_ts_monitor;

MRC_CLASS_DECLARE(mrc_ts, struct mrc_ts);

void mrc_ts_set_dt(struct mrc_ts *ts, float dt);
void mrc_ts_set_state(struct mrc_ts *ts, struct mrc_f1 *x);
void mrc_ts_add_monitor(struct mrc_ts *ts, struct mrc_ts_monitor *mon);
void mrc_ts_set_rhs_function(struct mrc_ts *ts,
			     void (*rhsf)(void *ctx, struct mrc_f1 *x,
					  struct mrc_f1 *rhs),
			     void *ctx);
void mrc_ts_solve(struct mrc_ts *ts);

#endif
