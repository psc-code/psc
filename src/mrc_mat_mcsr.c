
#include "mrc_mat_private.h"

#include <stdlib.h>

// ======================================================================
// mrc_mat "mcsr"

struct mrc_mat_mcsr_row {
  int idx;
  int first_entry;
};

struct mrc_mat_mcsr_entry {
  int idx;
  float val;
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

// ----------------------------------------------------------------------
// mrc_mat_mcsr_setup

static void
mrc_mat_mcsr_setup(struct mrc_mat *mat)
{
  struct mrc_mat_mcsr *sub = mrc_mat_mcsr(mat);

  sub->nr_rows_alloced = 1000;
  sub->nr_entries_alloced = 2000;

  sub->rows = calloc(sub->nr_rows_alloced, sizeof(*sub->rows));
  sub->entries = calloc(sub->nr_entries_alloced, sizeof(*sub->entries));

  sub->nr_entries = 0;
  sub->nr_rows = 0;
}

// ----------------------------------------------------------------------
// mrc_mat_mcsr_destroy

static void
mrc_mat_mcsr_destroy(struct mrc_mat *mat)
{
  struct mrc_mat_mcsr *sub = mrc_mat_mcsr(mat);

  free(sub->rows);
  free(sub->entries);
}

// ----------------------------------------------------------------------
// mrc_mat_mcsr_add_value

static void
mrc_mat_mcsr_add_value(struct mrc_mat *mat, int row_idx, int col_idx, float val)
{
  struct mrc_mat_mcsr *sub = mrc_mat_mcsr(mat);

  if (sub->nr_rows == 0 ||
      sub->rows[sub->nr_rows - 1].idx != row_idx) {
    // start new row
    if (sub->nr_rows >= sub->nr_rows_alloced - 1) {
      sub->nr_rows_alloced *= 2;
      sub->rows = realloc(sub->rows, sub->nr_rows_alloced * sizeof(*sub->rows));
    }
    sub->rows[sub->nr_rows].idx = row_idx;
    sub->rows[sub->nr_rows].first_entry = sub->nr_entries;
    sub->nr_rows++;
  }

  // if we already have an entry for this column in the current row, just add to it
  for (int i = sub->rows[sub->nr_rows - 1].first_entry; i < sub->nr_entries; i++) {
    if (sub->entries[i].idx == col_idx) {
      sub->entries[i].val += val;
      return;
    }
  }

  // otherwise, need to append a new entry
  if (sub->nr_entries >= sub->nr_entries_alloced) {
    sub->nr_entries_alloced *= 2;
    sub->entries = realloc(sub->entries, sub->nr_entries_alloced * sizeof(*sub->entries));
  }
  sub->entries[sub->nr_entries].idx = col_idx;
  sub->entries[sub->nr_entries].val = val;
  sub->nr_entries++;
}

// ----------------------------------------------------------------------
// mrc_mat_mcsr_assemble

static void
mrc_mat_mcsr_assemble(struct mrc_mat *mat)
{
  struct mrc_mat_mcsr *sub = mrc_mat_mcsr(mat);

  sub->rows[sub->nr_rows].first_entry = sub->nr_entries;
  mprintf("nr_rows %d nr_entries %d\n", sub->nr_rows, sub->nr_entries);
}

// ----------------------------------------------------------------------
// mrc_mat_mcsr_apply

static void
mrc_mat_mcsr_apply(struct mrc_fld *y, struct mrc_mat *mat, struct mrc_fld *x)
{
  struct mrc_mat_mcsr *sub = mrc_mat_mcsr(mat);

  assert(x->_size_of_type == sizeof(float));
  float *x_arr = x->_arr;
  float *y_arr = y->_arr;
    
  for (int row = 0; row < sub->nr_rows; row++) {
    int row_idx = sub->rows[row].idx;
    float sum = 0.;
    for (int entry = sub->rows[row].first_entry;
	 entry < sub->rows[row + 1].first_entry; entry++) {
      int col_idx = sub->entries[entry].idx;
      float val = sub->entries[entry].val;
      sum += val * x_arr[col_idx];
    }
    y_arr[row_idx] = sum;
  }
}

// ----------------------------------------------------------------------
// mrc_mat_mcsr_apply_in_place

static void
mrc_mat_mcsr_apply_in_place(struct mrc_mat *mat, struct mrc_fld *x)
{
  struct mrc_mat_mcsr *sub = mrc_mat_mcsr(mat);

  assert(x->_size_of_type == sizeof(float));
  float *arr = x->_arr;
    
  for (int row = 0; row < sub->nr_rows; row++) {
    int row_idx = sub->rows[row].idx;
    float sum = 0.;
    for (int entry = sub->rows[row].first_entry;
	 entry < sub->rows[row + 1].first_entry; entry++) {
      int col_idx = sub->entries[entry].idx;
      float val = sub->entries[entry].val;
      sum += val * arr[col_idx];
    }
    arr[row_idx] = sum;
  }
}

// ----------------------------------------------------------------------
// mrc_mat subclass "mcsr"

struct mrc_mat_ops mrc_mat_mcsr_ops = {
  .name                  = "mcsr",		
  .size                  = sizeof(struct mrc_mat_mcsr),
  .setup                 = mrc_mat_mcsr_setup,
  .destroy               = mrc_mat_mcsr_destroy,
  .add_value             = mrc_mat_mcsr_add_value,
  .assemble              = mrc_mat_mcsr_assemble,
  .apply                 = mrc_mat_mcsr_apply,
  .apply_in_place        = mrc_mat_mcsr_apply_in_place,
};

