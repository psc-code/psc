#include "mrc_mat_private.h"
#include "mrc_fld_as_double.h"
#include "mrc_ddc_private.h" 
#include <stdlib.h>

// ======================================================================
// mrc_mat "mcsr"

struct mrc_mat_mcsr_row {
  int idx;
  int first_entry;
};

struct mrc_mat_mcsr_entry {
  int idx;
  mrc_fld_data_t val;
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

  // mrc_mat "mcsr" works on a single proc only
  int size;
  MPI_Comm_size(mrc_mat_comm(mat), &size);
  assert(size == 1);

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
//
// WARNING, all elements for any given row must be added contiguously!

static void
mrc_mat_mcsr_realloc_rows_if_needed(struct mrc_mat *mat, int new_row)
{
  struct mrc_mat_mcsr *sub = mrc_mat_mcsr(mat);

  // need one final row for the final index into the entries array
  if (new_row >= sub->nr_rows_alloced - 1) {
    sub->nr_rows_alloced *= 2;
    sub->rows = realloc(sub->rows, sub->nr_rows_alloced * sizeof(*sub->rows));
  }
}

static void
mrc_mat_mcsr_add_value(struct mrc_mat *mat, int row_idx, int col_idx, double val)
{
  struct mrc_mat_mcsr *sub = mrc_mat_mcsr(mat);

  assert(row_idx >= 0 && row_idx < mat->m);
  assert(col_idx >= 0 && col_idx < mat->n);

  if (sub->nr_rows == 0 ||
      sub->rows[sub->nr_rows - 1].idx != row_idx) {
    // start new row
    mrc_mat_mcsr_realloc_rows_if_needed(mat, sub->nr_rows);
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
  mrc_mat_mcsr_realloc_rows_if_needed(mat, 0);
  sub->rows[sub->nr_rows].first_entry = sub->nr_entries;
}

// FIXME: semantics are different (wrong!) if empty rows are present:
// In mcsr, output vector values at empty rows are left unchanged, rather than zeroed

// ----------------------------------------------------------------------
// mrc_mat_mcsr_print

static void
mrc_mat_mcsr_print(struct mrc_mat *mat)
{
  struct mrc_mat_mcsr *sub = mrc_mat_mcsr(mat);

  //  mprintf("nr_rows = %d\n", sub->nr_rows);
  for (int row = 0; row < sub->nr_rows; row++) {
    int row_idx = sub->rows[row].idx;
    //    mprintf("row_idx = %d first_entry %d \n", row_idx, sub->rows[row].first_entry);
    for (int entry = sub->rows[row].first_entry;
	 entry < sub->rows[row + 1].first_entry; entry++) {
      int col_idx = sub->entries[entry].idx;
      mrc_fld_data_t val = sub->entries[entry].val;
      mprintf("row %d col %d val %g\n", row_idx, col_idx, val);
    }
  }
}

// ----------------------------------------------------------------------
// mrc_mat_mcsr_apply

static void
mrc_mat_mcsr_apply(struct mrc_fld *y, struct mrc_mat *mat, struct mrc_fld *x)
{
  struct mrc_mat_mcsr *sub = mrc_mat_mcsr(mat);

  assert(x->_size_of_type == sizeof(mrc_fld_data_t));
  mrc_fld_data_t *x_arr = x->_arr;
  mrc_fld_data_t *y_arr = y->_arr;
    
  for (int row = 0; row < sub->nr_rows; row++) {
    int row_idx = sub->rows[row].idx;
    mrc_fld_data_t sum = 0.;
    for (int entry = sub->rows[row].first_entry;
	 entry < sub->rows[row + 1].first_entry; entry++) {
      int col_idx = sub->entries[entry].idx;
      mrc_fld_data_t val = sub->entries[entry].val;
      sum += val * x_arr[col_idx];
    }
    y_arr[row_idx] = sum;
  }
}

// ----------------------------------------------------------------------
// mrc_mat_mcsr_apply_add

static void
mrc_mat_mcsr_apply_add(struct mrc_fld *y, struct mrc_mat *mat, struct mrc_fld *x)
{
  struct mrc_mat_mcsr *sub = mrc_mat_mcsr(mat);

  assert(x->_size_of_type == sizeof(mrc_fld_data_t));
  mrc_fld_data_t *x_arr = x->_arr;
  mrc_fld_data_t *y_arr = y->_arr;
    
  for (int row = 0; row < sub->nr_rows; row++) {
    int row_idx = sub->rows[row].idx;
    mrc_fld_data_t sum = 0.;
    for (int entry = sub->rows[row].first_entry;
	 entry < sub->rows[row + 1].first_entry; entry++) {
      int col_idx = sub->entries[entry].idx;
      mrc_fld_data_t val = sub->entries[entry].val;
      sum += val * x_arr[col_idx];
    }
    // FIXME, the only difference to apply() is the "+" in "+=", should be consolidated
    y_arr[row_idx] += sum;
  }
}

// ----------------------------------------------------------------------
// mrc_mat_mcsr_apply_in_place

static void
mrc_mat_mcsr_apply_in_place(struct mrc_mat *mat, struct mrc_fld *x)
{
  
  struct mrc_mat_mcsr *sub = mrc_mat_mcsr(mat);
  struct mrc_domain *domain = x->_domain; 
  struct mrc_ddc *ddc = mrc_domain_get_ddc(domain);

  int len = x->_len;
  if ( ddc->size == 1) {
      assert(x->_size_of_type == sizeof(mrc_fld_data_t));
      mrc_fld_data_t *arr = x->_arr;
      for (int row = 0; row < sub->nr_rows; row++) {
	int row_idx = sub->rows[row].idx;
	assert(row_idx < len);
	mrc_fld_data_t sum = 0.;
	for (int entry = sub->rows[row].first_entry;
	     entry < sub->rows[row + 1].first_entry; entry++) {
	  int col_idx = sub->entries[entry].idx;
	  assert(col_idx < len);
	  mrc_fld_data_t val = sub->entries[entry].val;
	  sum += val * arr[col_idx];
	}
	arr[row_idx] = sum;
      } 
 
  } else {    
    MPI_Comm comm = mrc_domain_comm(domain);   
    int rn[ddc->size]; 
    int ds[ddc->size]; 
    MPI_Allgather( &mat->n, 1, MPI_INT, rn, 1, MPI_INT, comm);
    int sz=0;
    for (int jj=0; jj < ddc->size; jj++) 
      {ds[jj]=sz; sz+=rn[jj];} 

    int nlpatches = mrc_fld_nr_patches(x); 
    int ngpatches; 
    mrc_domain_get_nr_global_patches(domain, &ngpatches); 
    int nppatch = mat->n/nlpatches;
    struct mrc_patch_info info; 
    // loop through all patches to find the first one 
    // assigned to this rank
    for (int pp=0; pp<ngpatches-1; pp++) { 
      mrc_domain_get_global_patch_info( domain, pp, &info); 
      if (info.rank == ddc->rank) {	
	break;
      }     
    }
    // allocate memory for gather vector 
    // each rank will recieve full x->_arr
    // this might not scale well but is easy to implement 
    mrc_fld_data_t *arr1= x->_arr; 
    mrc_fld_data_t *loc_x = (mrc_fld_data_t*) calloc(sz,sizeof(mrc_fld_data_t));   
    // FIXME: MPI_DOUBLE hardcoded here. 
    MPI_Allgatherv( &(arr1[0]), mat->n, MPI_DOUBLE, &(loc_x[0]), 
		    rn, ds, MPI_DOUBLE, comm);   
    mrc_fld_data_t *arr = loc_x;
    for (int row = 0; row < sub->nr_rows; row++) {
	int row_idx = sub->rows[row].idx;
	mrc_fld_data_t sum = 0.;
	for (int entry = sub->rows[row].first_entry;
	     entry < sub->rows[row + 1].first_entry; entry++) {
	  int col_idx = sub->entries[entry].idx;
	  mrc_fld_data_t val = sub->entries[entry].val;
	  sum += val * arr[col_idx];
	  } 
	if (ddc->rank !=0 ) {		  
	  arr1[row_idx - (nppatch*(info.global_patch))]=sum;
	} else {
	  arr1[row_idx]=sum;
	}
    }  
    free(loc_x);    
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
  .apply_add             = mrc_mat_mcsr_apply_add,
  .apply_in_place        = mrc_mat_mcsr_apply_in_place,
  .print                 = mrc_mat_mcsr_print,
};

