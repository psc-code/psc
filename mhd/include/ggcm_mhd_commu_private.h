
#ifndef GGCM_MHD_COMMU_PRIVATE_H
#define GGCM_MHD_COMMU_PRIVATE_H

#include "ggcm_mhd_commu.h"

struct ggcm_mhd_commu {
  struct mrc_obj obj;

  struct ggcm_mhd *mhd;
};

struct ggcm_mhd_commu_ops {
  MRC_SUBCLASS_OPS(struct ggcm_mhd_commu);
  void (*run)(struct ggcm_mhd_commu *commu, int mb, int me);
};

extern struct ggcm_mhd_commu_ops ggcm_mhd_commu_c_ops;

#define ggcm_mhd_commu_ops(commu) ((struct ggcm_mhd_commu_ops *)((commu)->obj.ops))

#endif
