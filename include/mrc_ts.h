
#ifndef MRC_TS_H
#define MRC_TS_H

#include <mrc_obj.h>

#include <stdio.h>

struct mrc_ts_monitor;

MRC_CLASS_DECLARE(mrc_ts, struct mrc_ts);

float mrc_ts_time(struct mrc_ts *ts);
float mrc_ts_dt(struct mrc_ts *ts);
int mrc_ts_step_number(struct mrc_ts *ts);
void mrc_ts_set_time(struct mrc_ts *ts, float time);
void mrc_ts_set_dt(struct mrc_ts *ts, float dt);
void mrc_ts_set_step_number(struct mrc_ts *ts, int n);
void mrc_ts_set_solution(struct mrc_ts *ts, struct mrc_obj *x);
void mrc_ts_add_monitor(struct mrc_ts *ts, struct mrc_ts_monitor *mon);
void mrc_ts_set_context(struct mrc_ts *ts, struct mrc_obj *ctx_obj);
void mrc_ts_set_get_dt_function(struct mrc_ts *ts,
				double (*get_dt_f)(void *ctx, struct mrc_ts *ts,
						   struct mrc_obj *x),
				void *ctx);
void mrc_ts_set_rhs_function(struct mrc_ts *ts,
			     void (*rhsf)(void *ctx, struct mrc_obj *x, float time,
					  struct mrc_obj *rhs),
			     void *ctx);
void mrc_ts_set_step_function(struct mrc_ts *ts,
			      void (*stepf)(void *ctx, struct mrc_ts *ts,
					    struct mrc_obj *x),
			      void *ctx);
void mrc_ts_set_pre_step_function(struct mrc_ts *ts,
				  void (*pre_step)(void *ctx, struct mrc_ts *ts,
						     struct mrc_obj *x),
				  void *ctx);
void mrc_ts_set_post_step_function(struct mrc_ts *ts,
				  void (*post_step)(void *ctx, struct mrc_ts *ts,
						    struct mrc_obj *x),
				  void *ctx);

void mrc_ts_step(struct mrc_ts *ts);
void mrc_ts_solve(struct mrc_ts *ts);

struct mrc_ts *mrc_ts_create_std(MPI_Comm comm,
				 void (*diagf)(void *ctx, float time,
					       struct mrc_obj *x, FILE *file),
				 void *diagf_ctx);

#endif
