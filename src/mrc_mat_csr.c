
#include "mrc_mat_private.h"

#include "mrc_fld_as_double.h" // FIXME, has to remain double, otherwise won't match mrc_mat_private.h
#include "mrc_vec.h"
#include "mrc_bits.h"

#include <stdlib.h>
#include <string.h>

// ----------------------------------------------------------------------
// mrc_mat_csr_create

static void
mrc_mat_csr_create(struct mrc_mat *mat)
{
  struct mrc_mat_csr *sub = mrc_mat_csr(mat);
  sub->vals = mrc_vec_create(mrc_mat_comm(mat));
  sub->cols = mrc_vec_create(mrc_mat_comm(mat));
  sub->rows = mrc_vec_create(mrc_mat_comm(mat));
  
  sub->_init_vals = NULL;
  sub->_init_cols = NULL;
  sub->_nr_cols = NULL;
  sub->_nr_cols_alloced = NULL;
  
  sub->is_assembled = false;
}

// ----------------------------------------------------------------------
// mrc_mat_csr_setup

static void
mrc_mat_csr_setup(struct mrc_mat *mat)
{
  struct mrc_mat_csr *sub = mrc_mat_csr(mat);

  // mrc_mat "csr" works on a single proc only
  int size;
  MPI_Comm_size(mrc_mat_comm(mat), &size);
  assert(size == 1);
  
  assert(sub->nr_initial_cols > 0);
  
  sub->nr_vals = 0;
  sub->nr_rows = mat->m;

  sub->_init_vals = calloc(sub->nr_rows, sizeof(*sub->_init_vals));
  sub->_init_cols = calloc(sub->nr_rows, sizeof(*sub->_init_vals));
  sub->_nr_cols = calloc(sub->nr_rows, sizeof(*sub->_nr_cols));
  sub->_nr_cols_alloced = calloc(sub->nr_rows, sizeof(*sub->_nr_cols_alloced));
  
  sub->_nr_rows_alloced = 0;
}

// ----------------------------------------------------------------------
// mrc_mat_csr_destroy

static void
mrc_mat_csr_destroy(struct mrc_mat *mat)
{
  struct mrc_mat_csr *sub = mrc_mat_csr(mat);
  
  if (mrc_mat_is_setup(mat) && !sub->is_assembled) {
    assert(sub->_init_vals != NULL && sub->_init_cols != NULL);
    for (int i=0; i < sub->nr_rows; i++) {
      if (sub->_init_vals[i] != NULL) {
        assert(sub->_init_cols[i] != NULL);
        free(sub->_init_vals[i]);
        free(sub->_init_cols[i]);
      } else {
        assert(sub->_init_cols[i] == NULL);
      }
    }
    free(sub->_init_vals);
    free(sub->_init_cols);
    
    assert(sub->_nr_cols != NULL && sub->_nr_cols_alloced != NULL);
    free(sub->_nr_cols);
    free(sub->_nr_cols_alloced);
    
    sub->_init_vals = NULL;
    sub->_init_cols = NULL;
    sub->_nr_cols = NULL;
    sub->_nr_cols_alloced = NULL;
  } else {
    // vals, cols, and rows are mrc_vecs and in _desr, so they're auto-cleaned
    assert(sub->_init_vals == NULL && sub->_init_cols == NULL &&
           sub->_nr_cols == NULL && sub->_nr_cols_alloced == NULL);
  }
}

// ----------------------------------------------------------------------
// _mrc_mat_csr_expand_row_if_needed

static void
_mrc_mat_csr_expand_row_if_needed(struct mrc_mat *mat, int row_idx)
{
  struct mrc_mat_csr *sub = mrc_mat_csr(mat);
  
  if (sub->_init_vals[row_idx] == NULL) {
    sub->_init_vals[row_idx] = calloc(sub->nr_initial_cols, sizeof(*sub->_init_vals));
    sub->_init_cols[row_idx] = calloc(sub->nr_initial_cols, sizeof(*sub->_init_cols));
    sub->_nr_vals_alloced += sub->nr_initial_cols;
    sub->_nr_cols_alloced[row_idx] = sub->nr_initial_cols;
    sub->_nr_rows_alloced += 1;
  } else if (sub->_nr_cols[row_idx] + 1 >= sub->_nr_cols_alloced[row_idx]) {
    // if there is not more room for another value in a given row, then
    // double the memory allocation for that row
    int new_size;
    int additional_cols = sub->growth_factor * sub->_nr_cols_alloced[row_idx];
    if (additional_cols == 0) {
      additional_cols = 1;
    }
    new_size = sub->_nr_cols_alloced[row_idx] + additional_cols;
    // mprintf("* expand_row realloc: row %d  %d -> %d\n",
    //         row_idx, sub->_nr_cols_alloced[row_idx], new_size);
    sub->_init_vals[row_idx] = realloc(sub->_init_vals[row_idx],
                                       new_size * sizeof(*sub->_init_vals[row_idx]));
    sub->_init_cols[row_idx] = realloc(sub->_init_cols[row_idx],
                                       new_size * sizeof(*sub->_init_cols[row_idx]));
    sub->_nr_vals_alloced += additional_cols;
    sub->_nr_cols_alloced[row_idx] = new_size;
  }
  // make sure my logic was correct
  assert(sub->_nr_cols[row_idx] + 1 <= sub->_nr_cols_alloced[row_idx]);
}

// ----------------------------------------------------------------------
// mrc_mat_csr_add_value

static void
mrc_mat_csr_add_value(struct mrc_mat *mat, int row_idx, int col_idx, double val)
{
  struct mrc_mat_csr *sub = mrc_mat_csr(mat);

  assert(row_idx >= 0 && row_idx < mat->m);
  assert(col_idx >= 0 && col_idx < mat->n);

  assert(!sub->is_assembled);
  
  int nr_cols = sub->_nr_cols[row_idx];

  // check if there's already a value in that location of the matrix
  if (sub->_init_vals[row_idx] != NULL) {
    for (int icol=0; icol < sub->_nr_cols[row_idx]; icol++) {
      if (sub->_init_cols[row_idx][icol] == col_idx) {
        // if the's already a value in this position, just add to it
        sub->_init_vals[row_idx][icol] += val;
        if (sub->_init_vals[row_idx][icol] == 0.0) {
          // the value has gone to 0, remove it from the matrix
          // mprintf("removing value: row %d  col %d  val %lg\n",
          //         row_idx, col_idx, sub->_init_vals[row_idx][icol]);
          for (int j=icol; j < nr_cols - 1; j++) {
            sub->_init_vals[row_idx][j] = sub->_init_vals[row_idx][j + 1];
            sub->_init_cols[row_idx][j] = sub->_init_cols[row_idx][j + 1];
          }
          sub->_nr_cols[row_idx] -= 1;
          sub->nr_vals -= 1;
        }
        return;
      } else if (sub->_init_cols[row_idx][icol] > col_idx) {
        // bail early since we know that value is currently 0
        break;
      }
    }
  }
  
  // just in case
  if (val == 0.0) {
    return;
  }
  _mrc_mat_csr_expand_row_if_needed(mat, row_idx);

  // if (sub->_nr_cols[row_idx] > 0 &&
  //     col_idx < sub->_init_cols[row_idx][sub->_nr_cols[row_idx] - 1]) {
  //   mprintf("WARNING: adding col out of order: (%d, %d) = %lg\n",
  //           row_idx, col_idx, val);
  // }

  // find insert index
  int i_insert = 0;
  for (int i=nr_cols; i > 0; i--){
    if (sub->_init_cols[row_idx][i - 1] < col_idx) {
      i_insert = i;
      break;
    } else {
      sub->_init_vals[row_idx][i] = sub->_init_vals[row_idx][i - 1];
      sub->_init_cols[row_idx][i] = sub->_init_cols[row_idx][i - 1];
    }
  }

  sub->_init_vals[row_idx][i_insert] = val;
  sub->_init_cols[row_idx][i_insert] = col_idx;
  sub->_nr_cols[row_idx] += 1;
  sub->nr_vals += 1;
}

// ----------------------------------------------------------------------
// mrc_mat_csr_assemble

static void
mrc_mat_csr_assemble(struct mrc_mat *mat)
{
  struct mrc_mat_csr *sub = mrc_mat_csr(mat);

  mrc_vec_set_type(sub->rows, "int");
  mrc_vec_set_param_int(sub->rows, "len", sub->nr_rows + 1);  
  mrc_vec_set_type(sub->vals, FLD_TYPE);
  mrc_vec_set_param_int(sub->vals, "len", sub->nr_vals);
  mrc_vec_set_type(sub->cols, "int");
  mrc_vec_set_param_int(sub->cols, "len", sub->nr_vals);
  mrc_vec_setup(sub->rows);
  mrc_vec_setup(sub->vals);
  mrc_vec_setup(sub->cols);

  mrc_fld_data_t *vals = mrc_vec_get_array(sub->vals);
  int *cols = mrc_vec_get_array(sub->cols);
  int *rows = mrc_vec_get_array(sub->rows);

  // if (sub->verbose) {
  // mprintf("!MAT_ASSEMBLE: %d -> %d  (-%d%%), %g MB -> %g MB\n",
  //         sub->_nr_vals_alloced, sub->nr_vals,
  //         (int)(100 * (1 - (1.0 * sub->nr_vals) / MAX(sub->_nr_vals_alloced, 1))),
  //         sizeof(mrc_fld_data_t) * sub->_nr_vals_alloced / 1e6,
  //         sizeof(mrc_fld_data_t) * sub->nr_vals / 1e6);
  // }
  if (sub->verbose && sub->nr_vals < 0.9 * sub->_nr_vals_alloced) {
    mprintf("NOTE: decreasing sparse matrix size by > 10%% on assemble: "
            "%d -> %d (-%d%%)\n", sub->_nr_vals_alloced, sub->nr_vals,
            (int)(100 * (1 - (1.0 * sub->nr_vals) / MAX(sub->_nr_vals_alloced, 1))));
  }

  int i = 0;
  for (int row=0; row < sub->nr_rows; row++) {
    int nr_cols = sub->_nr_cols[row];
    rows[row] = i;
    if (sub->_init_vals[row] != NULL) {
      assert(sub->_init_cols[row] != NULL);
      memcpy(&vals[i], sub->_init_vals[row], nr_cols * sizeof(*vals));
      memcpy(&cols[i], sub->_init_cols[row], nr_cols * sizeof(*cols));
      free(sub->_init_vals[row]);
      free(sub->_init_cols[row]);
      i += nr_cols;
    }
  }  
  rows[sub->nr_rows] = i;
  
  assert(sub->_nr_cols != NULL && sub->_nr_cols_alloced != NULL);
  free(sub->_init_vals);
  free(sub->_init_cols);
  free(sub->_nr_cols);
  free(sub->_nr_cols_alloced);
  sub->_init_vals = NULL;
  sub->_init_cols = NULL;
  sub->_nr_cols = NULL;
  sub->_nr_cols_alloced = NULL;
  
  sub->_nr_rows_alloced = 0;
  sub->_nr_vals_alloced = 0;
  
  mrc_vec_put_array(sub->vals, vals);
  mrc_vec_put_array(sub->cols, cols);
  mrc_vec_put_array(sub->rows, rows);
  sub->is_assembled = true;
}

// ----------------------------------------------------------------------
// mrc_mat_csr_print

static void
mrc_mat_csr_print(struct mrc_mat *mat)
{
  struct mrc_mat_csr *sub = mrc_mat_csr(mat);

  if (sub->is_assembled) {
    mrc_fld_data_t *vals = mrc_vec_get_array(sub->vals);
    int *cols = mrc_vec_get_array(sub->cols);
    int *rows = mrc_vec_get_array(sub->rows);
    for (int row_idx = 0; row_idx < sub->nr_rows; row_idx++) {
      for(int i=rows[row_idx]; i < rows[row_idx + 1]; i++) {
        mprintf("row %d col %d val %g\n", row_idx, cols[i], vals[i]);
      }
    }
    mrc_vec_put_array(sub->vals, vals);
    mrc_vec_put_array(sub->cols, cols);
    mrc_vec_put_array(sub->rows, rows);
  } else {
    for (int row_idx = 0; row_idx < sub->nr_rows; row_idx++) {
      for (int icol=0; icol < sub->_nr_cols[row_idx]; icol++) {
        mprintf("row %d col %d val %g\n", row_idx,
                sub->_init_cols[row_idx][icol], sub->_init_vals[row_idx][icol]);
      }
    }
  }
}

// ----------------------------------------------------------------------
// mrc_mat_csr_apply_gemv
// z = alpha * mat * x + beta * y

static inline void
_mrc_mat_csr_apply_gemv(struct mrc_vec *z, mrc_fld_data_t alpha,
                        struct mrc_mat *mat, struct mrc_vec *x,
                        mrc_fld_data_t beta, struct mrc_vec *y)
{
  struct mrc_mat_csr *sub = mrc_mat_csr(mat);
  const bool non_zero_beta = beta != 0.0;
  
  assert(mrc_vec_size_of_type(x) == sizeof(mrc_fld_data_t));
  assert(mrc_vec_size_of_type(y) == sizeof(mrc_fld_data_t));
  assert(mrc_vec_size_of_type(z) == sizeof(mrc_fld_data_t));
  assert(mrc_vec_len(x) == mat->n);
  assert(mrc_vec_len(y) == mat->m);
  assert(mrc_vec_len(z) == mrc_vec_len(y));

  assert(sub->is_assembled);

  mrc_fld_data_t *vals = mrc_vec_get_array(sub->vals);
  int *cols = mrc_vec_get_array(sub->cols);
  int *rows = mrc_vec_get_array(sub->rows);
  
  mrc_fld_data_t *x_arr = mrc_vec_get_array(x);
  mrc_fld_data_t *y_arr = mrc_vec_get_array(y);
  mrc_fld_data_t *z_arr = mrc_vec_get_array(z);

  for (int row_idx=0; row_idx < sub->nr_rows; row_idx++) {
    mrc_fld_data_t sum = 0.0;
    for (int i=rows[row_idx]; i < rows[row_idx + 1]; i++) {
      int col_idx = cols[i];
      sum += alpha * vals[i] * x_arr[col_idx];
    }
    if (non_zero_beta) {
      // this is protected by the if statement
      // in case y_arr == NAN but beta == 0.0
      sum += beta * y_arr[row_idx];
    }
    z_arr[row_idx] = sum;
  }

  mrc_vec_put_array(x, x_arr);
  mrc_vec_put_array(y, y_arr);
  mrc_vec_put_array(z, z_arr);
  mrc_vec_put_array(sub->vals, vals);
  mrc_vec_put_array(sub->cols, cols);
  mrc_vec_put_array(sub->rows, rows);  
}

// ----------------------------------------------------------------------
// mrc_mat_csr_apply
// y = mat * x

static void
mrc_mat_csr_apply(struct mrc_vec *y, struct mrc_mat *mat, struct mrc_vec *x)
{
  _mrc_mat_csr_apply_gemv(y, 1.0, mat, x, 0.0, y);
}

// ----------------------------------------------------------------------
// mrc_mat_csr_apply_in_place
// x = mat * x

static void
mrc_mat_csr_apply_in_place(struct mrc_mat *mat, struct mrc_vec *x)
{
  mrc_mat_csr_apply(x, mat, x);
}

// ----------------------------------------------------------------------
// mrc_mat_csr_apply_add
// y = mat * x + y

static void
mrc_mat_csr_apply_add(struct mrc_vec *y, struct mrc_mat *mat, struct mrc_vec *x)
{
  _mrc_mat_csr_apply_gemv(y, 1.0, mat, x, 1.0, y);
}

// ----------------------------------------------------------------------
// mrc_mat_csr_apply_general
// z = alpha * mat * x + beta * y

static void
mrc_mat_csr_apply_general(struct mrc_vec *z, double alpha,
                          struct mrc_mat *mat, struct mrc_vec *x,
                          double beta, struct mrc_vec *y)
{
  _mrc_mat_csr_apply_gemv(z, alpha, mat, x, beta, y);
}


// ----------------------------------------------------------------------
// mrc_mat_csr description
//
// nr_initial_* are just guesses so we don't have to do so much
// realloc'ing / moving data around at add_value and assemble

#define VAR(x) (void *)offsetof(struct mrc_mat_csr, x)
static struct param mrc_mat_csr_descr[] = {
  { "vals"             , VAR(vals)             , PARAM_OBJ(mrc_vec) },
  { "cols"             , VAR(cols)             , PARAM_OBJ(mrc_vec) },
  { "rows"             , VAR(rows)             , PARAM_OBJ(mrc_vec) },
  { "verbose"          , VAR(verbose)          , PARAM_BOOL(false)  },
  { "nr_initial_cols"  , VAR(nr_initial_cols)  , PARAM_INT(1),
  .help = "How many empty columnts to use for each new row" },
  { "growth_factor"    , VAR(growth_factor)    , PARAM_FLOAT(0.5),
  .help = "Fraction of current nr_cols to add to rows" },
  {},
};
#undef VAR

// ----------------------------------------------------------------------
// mrc_mat subclass "csr"

struct mrc_mat_ops mrc_mat_csr_ops = {
  .name                  = "csr",
  .size                  = sizeof(struct mrc_mat_csr),
  .param_descr           = mrc_mat_csr_descr,
  .create                = mrc_mat_csr_create,
  .setup                 = mrc_mat_csr_setup,
  .destroy               = mrc_mat_csr_destroy,
  .add_value             = mrc_mat_csr_add_value,
  .assemble              = mrc_mat_csr_assemble,
  .apply                 = mrc_mat_csr_apply,
  .apply_in_place        = mrc_mat_csr_apply_in_place,
  .apply_add             = mrc_mat_csr_apply_add,
  .apply_general         = mrc_mat_csr_apply_general,
  .print                 = mrc_mat_csr_print,
};
