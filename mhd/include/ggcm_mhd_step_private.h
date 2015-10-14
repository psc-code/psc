
#ifndef GGCM_MHD_STEP_PRIVATE_H
#define GGCM_MHD_STEP_PRIVATE_H

#include "ggcm_mhd_step.h"

#define MRC_FLD_CACHE_SIZE (20)
#define MRC_FLD_CACHE_COMPS (8)

#include "ggcm_mhd_diag_item.h"

struct mrc_fld_cache {
  int n;
  struct mrc_fld *flds[MRC_FLD_CACHE_SIZE];
};

struct ggcm_mhd_step {
  struct mrc_obj obj;

  struct ggcm_mhd *mhd;

  bool do_nwst; // calculate new dt next timestep?
  bool debug_dump;
  int profile_every; // print out profiling info every so many steps
  bool legacy_dt_handling; // handle timestep update as in legacy Fortran code

  struct mrc_fld_cache cache_1d[MRC_FLD_CACHE_COMPS];
};

struct ggcm_mhd_step_ops {
  MRC_SUBCLASS_OPS(struct ggcm_mhd_step);
  void (*newstep)(struct ggcm_mhd_step *, float *dtn);
  void (*pred)(struct ggcm_mhd_step *);
  void (*corr)(struct ggcm_mhd_step *);
  void (*calc_rhs)(struct ggcm_mhd_step *step, struct mrc_fld *rhs,
		   struct mrc_fld *x);
  void (*get_e_ec)(struct ggcm_mhd_step *step, struct mrc_fld *E,
                   struct mrc_fld *x);
  void (*run)(struct ggcm_mhd_step *step, struct mrc_fld *x);
  void (*setup_flds)(struct ggcm_mhd_step *step);
  void (*diag_item_zmask_run)(struct ggcm_mhd_step *step,
			      struct ggcm_mhd_diag_item *item,
			      struct mrc_io *io, struct mrc_fld *f,
			      int diag_type, float plane);
  void (*diag_item_rmask_run)(struct ggcm_mhd_step *step,
			      struct ggcm_mhd_diag_item *item,
			      struct mrc_io *io, struct mrc_fld *f,
			      int diag_type, float plane);

  int task_pred_nl1;
  int task_corr_nl1;
  int task_corr1_nl1;
  int task_pred_const;
  int task_corr_const;
  int task_corr1_const;
};

void ggcm_mhd_step_run_predcorr(struct ggcm_mhd_step *step, struct mrc_fld *x);

#define ggcm_mhd_step_ops(step) ((struct ggcm_mhd_step_ops *)(step)->obj.ops)

extern struct ggcm_mhd_step_ops ggcm_mhd_step_cweno_ops;
extern struct ggcm_mhd_step_ops ggcm_mhd_step_c_float_ops;
extern struct ggcm_mhd_step_ops ggcm_mhd_step_c_double_ops;
extern struct ggcm_mhd_step_ops ggcm_mhd_step_c2_float_ops;
extern struct ggcm_mhd_step_ops ggcm_mhd_step_c2_double_ops;
extern struct ggcm_mhd_step_ops ggcm_mhd_step_c3_float_ops;
extern struct ggcm_mhd_step_ops ggcm_mhd_step_c3_double_ops;
extern struct ggcm_mhd_step_ops ggcm_mhd_step_vlct_ops;
extern struct ggcm_mhd_step_ops ggcm_mhd_step_vl_ops;

struct mrc_fld *ggcm_mhd_step_get_1d_fld(struct ggcm_mhd_step *step, int nr_comps);
void ggcm_mhd_step_put_1d_fld(struct ggcm_mhd_step *step, struct mrc_fld *f);

#endif
