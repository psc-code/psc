
#ifndef MRC_TS_MONITOR_H
#define MRC_TS_MONITOR_H

#include <mrc_obj.h>

#include <mrc_ts.h>

MRC_CLASS_DECLARE(mrc_ts_monitor, struct mrc_ts_monitor);

void mrc_ts_monitor_run(struct mrc_ts_monitor *mon, struct mrc_ts *ts);

void mrc_ts_monitor_diag_set_function(struct mrc_ts_monitor *mon,
				      void (*diagf)(void *ctx, float time, struct mrc_obj *x,
						    FILE *file),
				      void *diagf_ctx);

#endif
