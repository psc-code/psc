
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

// ======================================================================
// mrc_mat "mcsr"
//
// This should be private to mrc_mat_mcsr.c, but mrc_mat_mcsr_mpi.c breaks
// proper separation of layers and needs to look into it.

struct mrc_mat_mcsr_row {
  int idx;
  int first_entry;
};

struct mrc_mat_mcsr_entry {
  int idx;
  double val; // FIXME should be mrc_fld_data_t
};

struct mrc_mat_mcsr {
  struct mrc_mat_mcsr_row *rows;
  struct mrc_mat_mcsr_entry *entries;
  int nr_rows;
  int nr_entries;
  int nr_rows_alloced;
  int nr_entries_alloced;
};

#define mrc_mat_mcsr(mat) mrc_to_subobj(mat, struct mrc_mat_mcsr)

#endif
