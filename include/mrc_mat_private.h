
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
  void (*apply)(struct mrc_vec *y, struct mrc_mat *mat, struct mrc_vec *x);
  void (*apply_add)(struct mrc_vec *y, struct mrc_mat *mat, struct mrc_vec *x);
  void (*apply_in_place)(struct mrc_mat *mat, struct mrc_vec *x);
  void (*apply_general)(struct mrc_vec *z, double alpha,
                        struct mrc_mat *mat, struct mrc_vec *x,
                        double beta, struct mrc_vec *y);
  void (*print)(struct mrc_mat *mat);
};

extern struct mrc_mat_ops mrc_mat_csr_ops;
extern struct mrc_mat_ops mrc_mat_csr_mpi_ops;
extern struct mrc_mat_ops mrc_mat_csr_slow_ops;
extern struct mrc_mat_ops mrc_mat_csr_slow_mpi_ops;
extern struct mrc_mat_ops mrc_mat_mcsr_ops;
extern struct mrc_mat_ops mrc_mat_mcsr_mpi_ops;
extern struct mrc_mat_ops mrc_mat_petsc_ops;

// ======================================================================
// mrc_mat "csr"
//
// This should be private to mrc_mat_csr.c, but mrc_mat_csr_mpi.c breaks
// proper separation of layers and needs to look into it.

struct mrc_mat_csr {
  bool is_assembled;
  int nr_rows;
  int nr_vals;
  
  // these are only valid after calling assemble
  struct mrc_vec *vals;  // type == double, length == nr_vals
  struct mrc_vec *cols;  // type == int, length == nr_vals
  struct mrc_vec *rows;  // type == int, length == nr_rows + 1
  // double *vals;  // length == _nr_vals_alloced (== nr_vals after assemble)
  // int *cols;  // length == _nr_vals_alloced (== nr_vals after assemble)
  // int *rows;  // length == nr_rows + 1
  
  // these are only valid before calling assemble
  double **_init_vals; // length == nr_rows
  int **_init_cols; // length == nr_rows
  int *_nr_cols; // length == nr_rows
  int *_nr_cols_alloced; // length == nr_rows
  int _nr_rows_alloced;  // probably used for asserts
  int _nr_vals_alloced;  // probably used for asserts

  int nr_initial_cols;  // how much initial space to allocate for each new row
};

#define mrc_mat_csr(mat) mrc_to_subobj(mat, struct mrc_mat_csr)

// ======================================================================
// mrc_mat "csr_slow"
//
// This should be private to mrc_mat_csr_slow.c, but mrc_mat_csr_slow_mpi.c breaks
// proper separation of layers and needs to look into it.

struct mrc_mat_csr_slow {
  // in principle, these could be mrc_vecs, but there's a lot of
  // resizing that happens while the matrix is being filled
  double *vals;  // length == _nr_vals_alloced (== nr_vals after assemble)
  int *cols;  // length == _nr_vals_alloced (== nr_vals after assemble)
  int *rows;  // length == nr_rows + 1
  int nr_rows;
  int nr_vals;
  bool is_assembled;
  // nr_initial_* are just guesses so we don't have to do so much
  // realloc'ing / moving data around at add_value and assemble
  int nr_initial_vals;  // how much memory to allocate at creation time
  int nr_initial_cols;  // how much initial space to allocate for each new row  
  // _nr_cols is used while adding values since rows will have more space than
  // they need until assemble is called
  int *_nr_cols;  // length == nr_rows, after assemble this == rows[i+1]-rows[i]
  int _nr_vals_alloced;  // after assemble this == nr_vals
  // _nr_current_rows is used to minimize the overhead of adding new rows
  // if they are added in order. Note that before the matrix is assembled,
  // only rows < _nr_current_rows are valid
  int _nr_current_rows;  // after assemble this == nr_rows
};

#define mrc_mat_csr_slow(mat) mrc_to_subobj(mat, struct mrc_mat_csr_slow)

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
