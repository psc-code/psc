
#ifndef GGCM_MHD_STEP_PRIVATE_H
#define GGCM_MHD_STEP_PRIVATE_H

#include "ggcm_mhd_step.h"

struct ggcm_mhd_step {
  struct mrc_obj obj;

  struct ggcm_mhd *mhd;

  bool do_nwst; // calculate new dt next timestep?
  int profile_every; // print out profiling info every so many steps
};

struct ggcm_mhd_step_ops {
  MRC_SUBCLASS_OPS(struct ggcm_mhd_step);
  void (*pred)(struct ggcm_mhd_step *);
  void (*corr)(struct ggcm_mhd_step *);
  void (*calc_rhs)(struct ggcm_mhd_step *step, struct mrc_fld *rhs,
		   struct mrc_fld *x);
  void (*run)(struct ggcm_mhd_step *step, struct mrc_fld *x);

  int task_pred_nl1;
  int task_corr_nl1;
  int task_corr1_nl1;
  int task_pred_const;
  int task_corr_const;
  int task_corr1_const;

  int mhd_type; // works on fully vs semi-conservative state vector?
  const char *fld_type; // works on this kind of mrc_fld
  int nr_ghosts; // this is how many ghost points we need
};

void ggcm_mhd_step_run_predcorr(struct ggcm_mhd_step *step, struct mrc_fld *x);

#define ggcm_mhd_step_ops(step) ((struct ggcm_mhd_step_ops *)(step)->obj.ops)

extern struct ggcm_mhd_step_ops ggcm_mhd_step_cweno_ops;
extern struct ggcm_mhd_step_ops ggcm_mhd_step_c_float_ops;
extern struct ggcm_mhd_step_ops ggcm_mhd_step_c2_float_ops;


#endif
