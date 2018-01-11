
#ifndef MRC_TS_PRIVATE_H
#define MRC_TS_PRIVATE_H

#include <mrc_ts.h>

#include <stdio.h>

struct mrc_ts {
  struct mrc_obj obj;
  // parameters
  int max_steps;
  float max_time;
  double norm_time;
  double norm_time_scale;

  int n; // current timestep number
  float time; // current integration time
  float dt; // current dt
  double tnorm; // normalization factor to get I/O units
  struct mrc_obj *x; // current state vector
  struct mrc_obj *ctx_obj;
  void *get_dt_f_ctx;
  double (*get_dt_f)(void *ctx, struct mrc_ts *ts, struct mrc_obj *x);
  void *rhsf_ctx;
  void (*rhsf)(void *ctx, struct mrc_obj *rhs, float time, struct mrc_obj *x);
  struct mrc_obj *stepf_ctx;
  void (*stepf)(void *ctx, struct mrc_ts *ts, struct mrc_obj *x);
  struct mrc_obj *pre_step_ctx;
  void (*pre_step)(void *ctx, struct mrc_ts *ts, struct mrc_obj *x);
  struct mrc_obj *post_step_ctx;
  void (*post_step)(void *ctx, struct mrc_ts *ts, struct mrc_obj *x);

  struct mrc_obj *(*vec_duplicate)(struct mrc_obj *);
  void (*vec_copy)(struct mrc_obj *, struct mrc_obj *);
  void (*vec_axpy)(struct mrc_obj *, float, struct mrc_obj *);
  void (*vec_waxpy)(struct mrc_obj *, float, struct mrc_obj *, struct mrc_obj *);
  float (*vec_norm)(struct mrc_obj *);
  void (*vec_set)(struct mrc_obj *, float);

  list_t monitors;

  int nr_rhsf_evals; // statistics
};

// ----------------------------------------------------------------------
// helpers for subclasses

void mrc_ts_rhsf(struct mrc_ts *ts, struct mrc_obj *rhs, float time,
		 struct mrc_obj *x);
void mrc_ts_monitors_run(struct mrc_ts *ts);

// ----------------------------------------------------------------------
// basic vector ops

// FIXME, should be standard mrc_obj functionality
static inline struct mrc_obj *
mrc_ts_vec_duplicate(struct mrc_ts *ts, struct mrc_obj *vec)
{
  return ts->vec_duplicate(vec);
}

static inline void
mrc_ts_vec_copy(struct mrc_ts *ts, struct mrc_obj *vecx, struct mrc_obj *vecy)
{
  ts->vec_copy(vecx, vecy);
}

static inline void
mrc_ts_vec_axpy(struct mrc_ts *ts, struct mrc_obj *vecy, float alpha,
		struct mrc_obj *vecx)
{
  ts->vec_axpy(vecy, alpha, vecx);
}

static inline void
mrc_ts_vec_waxpy(struct mrc_ts *ts, struct mrc_obj *vecw, float alpha,
		 struct mrc_obj *vecx, struct mrc_obj *vecy)
{
  ts->vec_waxpy(vecw, alpha, vecx, vecy);
}

static inline float
mrc_ts_vec_norm(struct mrc_ts *ts, struct mrc_obj *vec)
{
  return ts->vec_norm(vec);
}

static inline void
mrc_ts_vec_set(struct mrc_ts *ts, struct mrc_obj *vec, float val)
{
  ts->vec_set(vec, val);
}

// ----------------------------------------------------------------------

struct mrc_ts_ops {
  MRC_SUBCLASS_OPS(struct mrc_ts);
  void (*step)(struct mrc_ts *);
  void (*solve)(struct mrc_ts *);
};


#ifdef HAVE_PETSC
#include <petscts.h>
// mrc_ts_petsc - struct layout needed by any monitor that hooks into petsc
// ----------------------------------------------------------------------

// libmrc doesn't understand the distinction between global and local
// vectors, but petsc does, and any implicit solver will require
// a global vector to come in, so we actually need to copy the state
// vector in a hacky localtoglobal and feed that in.
struct mrc_ts_petsc {
  TS petsc_ts; // Petsc timestepper object 
  struct mrc_fld *F; // Computed RHS function will be stored here
  Vec F_vec;
  struct mrc_fld *xg; // Global (no ghost) version of the state vector
  Vec xg_vec;
  struct mrc_fld *xlocal; // temp X local for use in the calc_rhs wrapper.
  Mat J; // Jacobian Matrix.
  Mat Pre; // Optional matrix from which to calculate the preconditioner
  int debug_rhs;
  int (*calc_jac)(Mat J, Vec x, float t, void *ctx);
  int (*get_jac_mat)(void *ctx, Mat *M);
  int (*calc_pre_mat)(Mat J, Vec x, float t, void *ctx);
  int (*get_pre_mat)(void *ctx, Mat *M);
  void *jacf_ctx;
  
  // whether or not the preconditioner is calculated from a different matrix
  bool sep_pre;

  // pointer to a petsc MatStructure flag (allocated in *setup*) which gives the structure
  // of the preconditioner matrix compared to the jacobian matrix. The value of the pointer targer 
  // **must* be set by the user when the structure of the preconditioner matric structure is known.
  // FIXME: we really ought to the give the calc jacobian functions access to the mrc_ts...
  // options are DIFFERENT_NONZERO_PATTERN, SUBSET_NONZERO_PATTERN, SAME_NONZERO_PATTERN
  void *pre_mat_structure;
};

extern struct mrc_ts_ops mrc_ts_petsc_ops;

#endif

extern struct mrc_ts_ops mrc_ts_step_ops;
extern struct mrc_ts_ops mrc_ts_rk2_ops;
extern struct mrc_ts_ops mrc_ts_rk4_ops;
extern struct mrc_ts_ops mrc_ts_rkf45_ops;
extern struct mrc_ts_ops mrc_ts_ode45_ops;


#endif
