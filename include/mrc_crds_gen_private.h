
#ifndef MRC_CRDS_GEN_PRIVATE_H
#define MRC_CRDS_GEN_PRIVATE_H

#include "mrc_crds_gen.h"

struct mrc_crds_gen {
  struct mrc_obj obj;
  // parameters
  struct mrc_crds *crds;
  int d;
  // To make multi-block domains works crds_gen can't rely on the
  // domain global dims
  int dims[3];
  
  int n;  // n interior points
  int sw;  // stencil width
  float xl;  // low edge
  float xh;  // high edge
};

struct mrc_crds_gen_ops {
  MRC_SUBCLASS_OPS(struct mrc_crds_gen);
  void (*run)(struct mrc_crds_gen *gen, double *xx, double *dx);
};

#define mrc_crds_gen_ops(gen) ((struct mrc_crds_gen_ops *)(gen)->obj.ops)

extern struct mrc_crds_gen_ops mrc_crds_gen_uniform_ops;
extern struct mrc_crds_gen_ops mrc_crds_gen_ggcm_x_tanh_ops;
extern struct mrc_crds_gen_ops mrc_crds_gen_ggcm_x_cubic_ops;
extern struct mrc_crds_gen_ops mrc_crds_gen_ggcm_yz_ops;
extern struct mrc_crds_gen_ops mrc_crds_gen_gaussian_ops;
extern struct mrc_crds_gen_ops mrc_crds_gen_two_gaussian_ops;

#endif

