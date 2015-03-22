
#include "mrc_mat_private.h"

#include "mrc_fld_as_double.h" // FIXME, has to remain double, otherwise won't match mrc_mat_private.h
#include "mrc_vec.h"

#include <stdlib.h>
#include <string.h>

// ----------------------------------------------------------------------
// mrc_mat_csr_slow_setup

static void
mrc_mat_csr_slow_setup(struct mrc_mat *mat)
{
  struct mrc_mat_csr_slow *sub = mrc_mat_csr_slow(mat);

  // mrc_mat "csr_slow" works on a single proc only
  int size;
  MPI_Comm_size(mrc_mat_comm(mat), &size);
  assert(size == 1);

  assert(sub->nr_initial_vals > 0);
  assert(sub->nr_initial_cols > 0);
  
  sub->nr_rows = mat->m;
  sub->_nr_current_rows = 0;
  sub->nr_vals = 0;
  sub->_nr_vals_alloced = sub->nr_initial_vals;
  sub->is_assembled = false;
  
  // size will not change
  sub->rows = calloc(sub->nr_rows + 1, sizeof(*sub->rows));
  sub->_nr_cols = calloc(sub->nr_rows, sizeof(*sub->_nr_cols));
  // size will probably change
  sub->vals = calloc(sub->_nr_vals_alloced, sizeof(*sub->vals));
  sub->cols = calloc(sub->_nr_vals_alloced, sizeof(*sub->cols));
}

// ----------------------------------------------------------------------
// mrc_mat_csr_slow_destroy

static void
mrc_mat_csr_slow_destroy(struct mrc_mat *mat)
{
  struct mrc_mat_csr_slow *sub = mrc_mat_csr_slow(mat);
  
  free(sub->rows);
  free(sub->_nr_cols);
  free(sub->vals);
  free(sub->cols);
}

// ----------------------------------------------------------------------
// mrc_mat_csr_slow_add_empty_rows_if_needed
//
// Fill the rows between sub->_nr_current_rows and nr_rows with empty rows.
// After this call, sub->rows up to and including nr_rows will be valid.
// Also, sub->_nr_cols up to and excluding nr_rows will be valid.

static void
_mrc_mat_csr_slow_expand_nr_rows_if_needed(struct mrc_mat *mat, int nr_rows)
{
  struct mrc_mat_csr_slow *sub = mrc_mat_csr_slow(mat);
  assert(nr_rows <= mat->m);
  
  if (nr_rows <= sub->_nr_current_rows) {
    return;
  }
  
  int i_last_row = sub->rows[sub->_nr_current_rows];
  for (int i=sub->_nr_current_rows; i < nr_rows; i++) {
    sub->rows[i] = i_last_row;
    sub->_nr_cols[i] = 0;
  }
  sub->rows[nr_rows] = i_last_row;
  sub->_nr_current_rows = nr_rows;
}

// ----------------------------------------------------------------------
// _mrc_mat_csr_slow_expand_row_if_needed
//
// If there isn't enough room for one more value in a row, double the
// number of columns for that row and move rows > row_idx accordingly
// If there isn't enough allocated memory to do that, then first
// allocate some more memory.

static void
_mrc_mat_csr_slow_expand_row_if_needed(struct mrc_mat *mat, int row_idx)
{
  struct mrc_mat_csr_slow *sub = mrc_mat_csr_slow(mat);

  int nr_cols_alloced = sub->rows[row_idx + 1] - sub->rows[row_idx];
  int nr_empty_cols = nr_cols_alloced - sub->_nr_cols[row_idx];
  assert(nr_empty_cols >= 0);
  // if we don't have space for any new values, grow the number of columns
  // in the current row by a factor of 2, then move the data of all subsequent
  // rows
  if (nr_empty_cols == 0) {
    int nr_new_cols;
    if (nr_cols_alloced == 0) {
      nr_new_cols = sub->nr_initial_cols;
    } else {
      nr_new_cols = nr_cols_alloced;
    }
    
    // do we need more memory to fit the new row size? if so, grow the size
    // of the allocation by a factor of 2
    int nr_vals_in_use = sub->rows[sub->_nr_current_rows];
    int nr_available_vals = sub->_nr_vals_alloced - nr_vals_in_use;
    if (nr_available_vals < nr_new_cols) {
      // NOTE: a factor of 2 seems too big for large matrices...
      int new_nr_vals_alloced = 2 * sub->_nr_vals_alloced;
      // new_nr_vals_alloced = sub->_nr_vals_alloced + sub->nr_initial_vals;

      // make sure the new allocation will actually fit the new cols
      int expansion_size = new_nr_vals_alloced - sub->_nr_vals_alloced;
      if (expansion_size < (nr_new_cols - nr_available_vals)) {
        expansion_size = nr_new_cols - nr_available_vals;
        new_nr_vals_alloced = sub->_nr_vals_alloced + expansion_size;
      }
      
      mprintf("Expanding matrix allocation: %d -> %d  (%+d)\n",
              sub->_nr_vals_alloced, new_nr_vals_alloced, expansion_size);
      sub->vals = realloc(sub->vals, new_nr_vals_alloced * sizeof(*sub->vals));
      sub->cols = realloc(sub->cols, new_nr_vals_alloced * sizeof(*sub->cols));
      sub->_nr_vals_alloced = new_nr_vals_alloced;
    }
    
    // now move the rows from row_idx to the end
    // mprintf("Expanding row %d: %d -> %d (%+d) and moving %d rows\n", row_idx,
    //         nr_cols_alloced, nr_cols_alloced + nr_new_cols, nr_new_cols,
    //         sub->_nr_current_rows - 1 - row_idx);
    sub->rows[sub->_nr_current_rows] += nr_new_cols;
    for (int i=sub->_nr_current_rows - 1; i > row_idx; i--) {
      // mprintf("moving row: %d\n", i);
      int n = sub->_nr_cols[i];
      int i_row = sub->rows[i];
      memmove(&sub->vals[i_row + nr_new_cols], &sub->vals[i_row], n * sizeof(*sub->vals));
      memmove(&sub->cols[i_row + nr_new_cols], &sub->cols[i_row], n * sizeof(*sub->cols));
      sub->rows[i] += nr_new_cols;
      // mprintf("row moved: %d\n", i);
    }
  }
}

// ----------------------------------------------------------------------
// mrc_mat_csr_slow_add_value

static void
mrc_mat_csr_slow_add_value(struct mrc_mat *mat, int row_idx, int col_idx, double val)
{
  struct mrc_mat_csr_slow *sub = mrc_mat_csr_slow(mat);

  // indices named *_idx are relative to the row/column (m/n) of the matrix
  // indices named i_* are relative to the value's position in sub->vals

  assert(row_idx >= 0 && row_idx < mat->m);
  assert(col_idx >= 0 && col_idx < mat->n);

  sub->is_assembled = false;
  
  if (row_idx < (sub->_nr_current_rows - 1)) {
    mprintf("WARNING:: adding rows out of order: (%d, %d) = %lg\n",
            row_idx, col_idx, val);
  }  
  _mrc_mat_csr_slow_expand_nr_rows_if_needed(mat, row_idx + 1);
  
  int nr_cols = sub->_nr_cols[row_idx];
  int i_next_col = sub->rows[row_idx] + nr_cols;

  if (col_idx < (nr_cols - 1)) {
    mprintf("WARNING:: adding cols out of order: (%d, %d) = %lg\n",
            row_idx, col_idx, val);
  }

  // check if there's already a value in that location of the matrix
  for (int i = sub->rows[row_idx]; i < i_next_col; i++) {
    if (col_idx == sub->cols[i]) {
      // if the's already a value in this position, just add to it
      sub->vals[i] += val;
      
      if (sub->vals[i] == 0.0) {
        // the value has gone to 0, remove it from the matrix
        for (int j=i; j < i_next_col - 1; j++){
          sub->vals[j] = sub->vals[j + 1];
          sub->cols[j] = sub->cols[j + 1];
        }
        assert(sub->_nr_cols[row_idx] > 0);  // this is probably not needed
        sub->_nr_cols[row_idx] -= 1;
        sub->nr_vals -= 1;
      }
      return;
    }
  }

  _mrc_mat_csr_slow_expand_row_if_needed(mat, row_idx);

  // int next_col_in_row_idx = sub->rows[row_idx] + sub->_nr_cols[row_idx];
  int i_insert = -1;

  if (nr_cols == 0){
    // this is the first column in the row
    i_insert = sub->rows[row_idx];
  } else {
    // find insert index
    i_insert = sub->rows[row_idx];
    for (int i=i_next_col; i > sub->rows[row_idx]; i--){
      if (sub->cols[i - 1] < col_idx) {
        i_insert = i;
        break;
      } else {
        sub->cols[i] = sub->cols[i - 1];
        sub->vals[i] = sub->vals[i - 1];
      }
    }
  }
  assert(i_insert >= 0);
  assert(i_insert < sub->rows[row_idx + 1]);
  // assert(i_insert < sub->nr_cols_alloced[row_idx]);  // should be done by expand_row_if_needed
  sub->vals[i_insert] = val;
  sub->cols[i_insert] = col_idx;
  sub->_nr_cols[row_idx] += 1;
  sub->nr_vals += 1;
}

// ----------------------------------------------------------------------
// mrc_mat_csr_slow_assemble

static void
mrc_mat_csr_slow_assemble(struct mrc_mat *mat)
{
  struct mrc_mat_csr_slow *sub = mrc_mat_csr_slow(mat);

  mprintf("starting matrix assemble\n");

  // compress the matrix so that there's no empty space at the end of rows
  for (int row_idx=1; row_idx < sub->_nr_current_rows + 1; row_idx++) {
    int end_prev_row = sub->rows[row_idx - 1] + sub->_nr_cols[row_idx - 1];
    int nr_empty_vals = sub->rows[row_idx] - end_prev_row;

    if (row_idx < sub->_nr_current_rows) {
      // move some data
      int i_next_col = sub->rows[row_idx] + sub->_nr_cols[row_idx];
      for (int i=sub->rows[row_idx]; i < i_next_col; i++) {
        sub->vals[i - nr_empty_vals] = sub->vals[i];
        sub->cols[i - nr_empty_vals] = sub->cols[i];
      }
    }
    sub->rows[row_idx] -= nr_empty_vals;
    assert(sub->rows[row_idx] ==
           sub->rows[row_idx - 1] + sub->_nr_cols[row_idx - 1]);
  }
  
  // now fill the matrix from _nr_nows_added to nr_rows with empty rows
  _mrc_mat_csr_slow_expand_nr_rows_if_needed(mat, sub->nr_rows);
  assert(sub->_nr_current_rows == sub->nr_rows);
  
  // mprintf("on assemble: nr_vals alloced %d -> %d  (%+d)\n",
  //         sub->_nr_vals_alloced, sub->rows[sub->nr_rows],
  //         sub->rows[sub->nr_rows] - sub->_nr_vals_alloced);
  sub->_nr_vals_alloced = sub->rows[sub->nr_rows];
  assert(sub->nr_vals == sub->_nr_vals_alloced);
  sub->vals = realloc(sub->vals, sub->nr_vals * sizeof(*sub->vals));
  sub->cols = realloc(sub->cols, sub->nr_vals * sizeof(*sub->cols));
  sub->is_assembled = true;
  mprintf("finished matrix assemble\n");
}

// ----------------------------------------------------------------------
// mrc_mat_csr_slow_print

static void
mrc_mat_csr_slow_print(struct mrc_mat *mat)
{
  struct mrc_mat_csr_slow *sub = mrc_mat_csr_slow(mat);

  //  mprintf("nr_rows = %d\n", sub->nr_rows);
  for (int row_idx = 0; row_idx < sub->nr_rows; row_idx++) {
    // using _nr_cols so the matrix doesn't have to be assembled yet
    int i_next_col = sub->rows[row_idx] + sub->_nr_cols[row_idx];    
    for (int i = sub->rows[row_idx]; i < i_next_col; i++){
      mprintf("row %d col %d val %g\n", row_idx, sub->cols[i], sub->vals[i]);
    }
  }
}

// ----------------------------------------------------------------------
// mrc_mat_csr_slow_apply_gemv
// z = alpha * mat * x + beta * y

static inline void
_mrc_mat_csr_slow_apply_gemv(struct mrc_vec *z, mrc_fld_data_t alpha,
                             struct mrc_mat *mat, struct mrc_vec *x,
                             mrc_fld_data_t beta, struct mrc_vec *y)
{
  struct mrc_mat_csr_slow *sub = mrc_mat_csr_slow(mat);
  const bool non_zero_alpha = alpha != 0.0;
  const bool non_zero_beta = beta != 0.0;
  
  assert(mrc_vec_size_of_type(x) == sizeof(mrc_fld_data_t));
  assert(mrc_vec_size_of_type(y) == sizeof(mrc_fld_data_t));
  assert(mrc_vec_size_of_type(z) == sizeof(mrc_fld_data_t));
  assert(mrc_vec_len(x) == mat->n);
  assert(mrc_vec_len(y) == mat->m);
  assert(mrc_vec_len(z) == mrc_vec_len(y));

  assert(sub->is_assembled);

  mrc_fld_data_t *x_arr = mrc_vec_get_array(x);
  mrc_fld_data_t *y_arr = mrc_vec_get_array(y);
  mrc_fld_data_t *z_arr = mrc_vec_get_array(z);

  for (int row_idx=0; row_idx < sub->nr_rows; row_idx++) {
    mrc_fld_data_t sum = 0.0;
    for (int i=sub->rows[row_idx]; i < sub->rows[row_idx + 1]; i++) {
      int col_idx = sub->cols[i];
      mrc_fld_data_t val = sub->vals[i];
      if (non_zero_alpha) {
        val *= alpha;
      }
      sum += val * x_arr[col_idx];
    }
    if (non_zero_beta) {
      sum += beta * y_arr[row_idx];
    }
    z_arr[row_idx] = sum;
  }

  mrc_vec_put_array(x, x_arr);
  mrc_vec_put_array(y, y_arr);
  mrc_vec_put_array(z, z_arr);
}


// ----------------------------------------------------------------------
// mrc_mat_csr_slow_apply
// y = mat * x

static void
mrc_mat_csr_slow_apply(struct mrc_vec *y, struct mrc_mat *mat, struct mrc_vec *x)
{
  _mrc_mat_csr_slow_apply_gemv(y, 1.0, mat, x, 0.0, y);
}

// ----------------------------------------------------------------------
// mrc_mat_csr_slow_apply_in_place
// x = mat * x

static void
mrc_mat_csr_slow_apply_in_place(struct mrc_mat *mat, struct mrc_vec *x)
{
  mrc_mat_csr_slow_apply(x, mat, x);
}

// ----------------------------------------------------------------------
// mrc_mat_csr_slow_apply_add
// y = mat * x + y

static void
mrc_mat_csr_slow_apply_add(struct mrc_vec *y, struct mrc_mat *mat, struct mrc_vec *x)
{
  _mrc_mat_csr_slow_apply_gemv(y, 1.0, mat, x, 1.0, y);
}

// ----------------------------------------------------------------------
// mrc_mat_csr_slow_apply_general
// z = alpha * mat * x + beta * y

static void
mrc_mat_csr_slow_apply_general(struct mrc_vec *z, double alpha,
                          struct mrc_mat *mat, struct mrc_vec *x,
                          double beta, struct mrc_vec *y)
{
  _mrc_mat_csr_slow_apply_gemv(z, alpha, mat, x, beta, y);
}


// ----------------------------------------------------------------------
// mrc_mat_csr_slow_mpi description
//
// nr_initial_* are just guesses so we don't have to do so much
// realloc'ing / moving data around at add_value and assemble

#define VAR(x) (void *)offsetof(struct mrc_mat_csr_slow, x)
static struct param mrc_mat_csr_slow_descr[] = {
  { "nr_initial_vals"        , VAR(nr_initial_vals)     , PARAM_INT(1),
  .help = "How many values to malloc at creation time" },
  { "nr_initial_cols"        , VAR(nr_initial_cols)     , PARAM_INT(1),
  .help = "How many empty columnts to use for each new row" },
  {},
};
#undef VAR

// ----------------------------------------------------------------------
// mrc_mat subclass "csr_slow"

struct mrc_mat_ops mrc_mat_csr_slow_ops = {
  .name                  = "csr_slow",
  .size                  = sizeof(struct mrc_mat_csr_slow),
  .param_descr           = mrc_mat_csr_slow_descr,
  .setup                 = mrc_mat_csr_slow_setup,
  .destroy               = mrc_mat_csr_slow_destroy,
  .add_value             = mrc_mat_csr_slow_add_value,
  .assemble              = mrc_mat_csr_slow_assemble,
  .apply                 = mrc_mat_csr_slow_apply,
  .apply_in_place        = mrc_mat_csr_slow_apply_in_place,
  .apply_add             = mrc_mat_csr_slow_apply_add,
  .apply_general         = mrc_mat_csr_slow_apply_general,
  .print                 = mrc_mat_csr_slow_print,
};
