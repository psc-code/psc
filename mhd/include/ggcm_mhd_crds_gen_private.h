
#ifndef GGCM_MHD_CRDS_GEN_PRIVATE_H
#define GGCM_MHD_CRDS_GEN_PRIVATE_H

#include "ggcm_mhd_crds_gen.h"

struct ggcm_mhd_crds_gen {
  struct mrc_obj obj;
  bool legacy_fd1;
};

struct ggcm_mhd_crds_gen_ops {
  MRC_SUBCLASS_OPS(struct ggcm_mhd_crds_gen);

  void (*run)(struct ggcm_mhd_crds_gen *gen, struct ggcm_mhd_crds *crds);
  void (*run_aux)(struct ggcm_mhd_crds_gen *gen, struct ggcm_mhd_crds *crds);
};

extern struct ggcm_mhd_crds_gen_ops ggcm_mhd_crds_gen_mrc_ops;

#endif
