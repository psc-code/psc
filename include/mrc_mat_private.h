
#ifndef MRC_MAT_PRIVATE_H
#define MRC_MAT_PRIVATE_H

#include "mrc_mat.h"

struct mrc_mat {
  struct mrc_obj obj;
};

struct mrc_mat_ops {
  MRC_SUBCLASS_OPS(struct mrc_mat);
  void (*add_value)(struct mrc_mat *mat, int row_idx, int col_idx, float val);
  void (*assemble)(struct mrc_mat *mat);
  void (*apply_in_place)(struct mrc_mat *mat, struct mrc_fld *x);
};

extern struct mrc_mat_ops mrc_mat_mcsr_ops;

#endif
