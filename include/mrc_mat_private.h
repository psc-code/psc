
#ifndef MRC_MAT_PRIVATE_H
#define MRC_MAT_PRIVATE_H

#include "mrc_mat.h"

struct mrc_mat {
  struct mrc_obj obj;

  // parameters
  int m; // number of local rows
  int n; // number of local columns
};

struct mrc_mat_ops {
  MRC_SUBCLASS_OPS(struct mrc_mat);
  void (*add_value)(struct mrc_mat *mat, int row_idx, int col_idx, double val);
  void (*assemble)(struct mrc_mat *mat);
  void (*apply)(struct mrc_fld *y, struct mrc_mat *mat, struct mrc_fld *x);
  void (*apply_add)(struct mrc_fld *y, struct mrc_mat *mat, struct mrc_fld *x);
  void (*apply_in_place)(struct mrc_mat *mat, struct mrc_fld *x);
  void (*print)(struct mrc_mat *mat);
};

extern struct mrc_mat_ops mrc_mat_mcsr_ops;
extern struct mrc_mat_ops mrc_mat_mcsr_mpi_ops;
extern struct mrc_mat_ops mrc_mat_petsc_ops;

#endif
