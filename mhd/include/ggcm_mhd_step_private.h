
#ifndef GGCM_MHD_STEP_PRIVATE_H
#define GGCM_MHD_STEP_PRIVATE_H

#include "ggcm_mhd_step.h"

struct ggcm_mhd_step {
  struct mrc_obj obj;

  struct ggcm_mhd *mhd;
};

struct ggcm_mhd_step_ops {
  MRC_SUBCLASS_OPS(struct ggcm_mhd_step);
  void (*pred)(struct ggcm_mhd_step *);
  void (*corr)(struct ggcm_mhd_step *);
  void (*calc_rhs)(struct ggcm_mhd_step *step, struct mrc_f3 *rhs,
		   struct mrc_f3 *x);

  int task_pred_nl1;
  int task_corr_nl1;
  int task_corr1_nl1;
  int task_pred_const;
  int task_corr_const;
  int task_corr1_const;
};

#define ggcm_mhd_step_ops(step) ((struct ggcm_mhd_step_ops *)(step)->obj.ops)

#endif
