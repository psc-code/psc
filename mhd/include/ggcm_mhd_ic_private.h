
#ifndef GGCM_MHD_IC_PRIVATE_H
#define GGCM_MHD_IC_PRIVATE_H

#include "ggcm_mhd_ic.h"

struct ggcm_mhd_ic {
  struct mrc_obj obj;
  struct ggcm_mhd *mhd;
};

struct ggcm_mhd_ic_ops {
  MRC_SUBCLASS_OPS(struct ggcm_mhd_ic);
  void (*run)(struct ggcm_mhd_ic *ic);
  double (*primitive)(struct ggcm_mhd_ic *ic, int m, double crd[3]);
  double (*primitive_bg)(struct ggcm_mhd_ic *ic, int m, double crd[3]);
  double (*vector_potential)(struct ggcm_mhd_ic *ic, int m, double crd[3]);
  double (*vector_potential_bg)(struct ggcm_mhd_ic *ic, int m, double crd[3]);
};

#define ggcm_mhd_ic_ops(ic) ((struct ggcm_mhd_ic_ops *)((ic)->obj.ops))

extern struct ggcm_mhd_ic_ops ggcm_mhd_ic_mirdip_float_ops;
extern struct ggcm_mhd_ic_ops ggcm_mhd_ic_mirdip_double_ops;
extern struct ggcm_mhd_ic_ops ggcm_mhd_ic_obstacle_double_ops;
extern struct ggcm_mhd_ic_ops ggcm_mhd_ic_gkeyll_ops;

#endif
