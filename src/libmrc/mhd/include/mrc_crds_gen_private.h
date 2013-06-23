
#ifndef MRC_CRDS_GEN_PRIVATE_H
#define MRC_CRDS_GEN_PRIVATE_H

#include "mrc_crds_gen.h"

struct mrc_crds_gen {
  struct mrc_obj obj;
  // parameters
  struct mrc_crds *crds;
  int d;
};

struct mrc_crds_gen_ops {
  MRC_SUBCLASS_OPS(struct mrc_crds_gen);
  void (*run)(struct mrc_crds_gen *gen, float *xx, float *dx);
};

#define mrc_crds_gen_ops(gen) ((struct mrc_crds_gen_ops *)(gen)->obj.ops)

extern struct mrc_crds_gen_ops mrc_crds_gen_ggcm_x_ops;
extern struct mrc_crds_gen_ops mrc_crds_gen_ggcm_yz_ops;

#endif

