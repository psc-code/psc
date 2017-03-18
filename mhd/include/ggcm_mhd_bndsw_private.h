
#ifndef GGCM_MHD_BNDSW_PRIVATE_H
#define GGCM_MHD_BNDSW_PRIVATE_H

#include "ggcm_mhd_bndsw.h"

struct ggcm_mhd_bndsw {
  struct mrc_obj obj;
  // parameters
  float alphafak;
};

struct ggcm_mhd_bndsw_ops {
  MRC_SUBCLASS_OPS(struct ggcm_mhd_bndsw);
  void (*at)(struct ggcm_mhd_bndsw *bndsw, float bntim, float xx[3],
	     float vals[]);
  void (*get_initial)(struct ggcm_mhd_bndsw *bndsw, float vals[]);
  void (*new_step)(struct ggcm_mhd_bndsw *bndsw);
};

extern struct ggcm_mhd_bndsw_ops ggcm_mhd_bndsw_none_ops;

#define ggcm_mhd_bndsw_ops(bndsw) ((struct ggcm_mhd_bndsw_ops *)((bndsw)->obj.ops))

#endif
