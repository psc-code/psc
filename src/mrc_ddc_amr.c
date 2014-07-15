
#include "mrc_ddc_private.h"

#include <mrc_domain.h>
#include <mrc_params.h>
#include <assert.h>
#include <stdlib.h>
#include <stdio.h>

// ======================================================================
// mrc_ddc_amr

struct mrc_ddc_amr_row {
  int patch;
  int idx;
  int first_entry;
};

struct mrc_ddc_amr_entry {
  int pidx;
  float val;
};

struct mrc_mat_mcsr {
  struct mrc_ddc_amr_row *rows;
  struct mrc_ddc_amr_entry *entries;
  int nr_rows;
  int nr_entries;
  int nr_rows_alloced;
  int nr_entries_alloced;
};

struct mrc_ddc_amr {
  struct mrc_mat_mcsr mat;

  struct mrc_domain *domain;
  int sw[3];
  int ib[3], im[4];
};

#define mrc_ddc_amr(ddc) mrc_to_subobj(ddc, struct mrc_ddc_amr)

// ----------------------------------------------------------------------
// mrc_ddc_amr_set_domain

static void
mrc_ddc_amr_set_domain(struct mrc_ddc *ddc, struct mrc_domain *domain)
{
  struct mrc_ddc_amr *sub = mrc_ddc_amr(ddc);
  sub->domain = domain;
}

// ----------------------------------------------------------------------
// mrc_ddc_amr_get_domain

static struct mrc_domain *
mrc_ddc_amr_get_domain(struct mrc_ddc *ddc)
{
  struct mrc_ddc_amr *sub = mrc_ddc_amr(ddc);
  return sub->domain;
}

// ----------------------------------------------------------------------
// mrc_ddc_amr_setup

static void
mrc_ddc_amr_setup(struct mrc_ddc *ddc)
{
  struct mrc_ddc_amr *sub = mrc_ddc_amr(ddc);
  struct mrc_mat_mcsr *mcsr = &sub->mat;
  assert(sub->domain);

  int ldims[3];
  mrc_domain_get_param_int3(sub->domain, "m", ldims);
  // needs to be compatible with how mrc_fld indexes its fields
  for (int d = 0; d < 3; d++) {
    sub->ib[d] = -sub->sw[d];
    sub->im[d] = ldims[d] + 2 * sub->sw[d];
  }

  mcsr->nr_rows_alloced = 1000;
  mcsr->nr_entries_alloced = 2000;

  mcsr->rows = calloc(mcsr->nr_rows_alloced, sizeof(*mcsr->rows));
  mcsr->entries = calloc(mcsr->nr_entries_alloced, sizeof(*mcsr->entries));

  mcsr->nr_entries = 0;
  mcsr->nr_rows = 0;
}

// ----------------------------------------------------------------------
// mrc_ddc_amr_destroy

static void
mrc_ddc_amr_destroy(struct mrc_ddc *ddc)
{
  struct mrc_ddc_amr *sub = mrc_ddc_amr(ddc);
  struct mrc_mat_mcsr *mcsr = &sub->mat;

  free(mcsr->rows);
  free(mcsr->entries);
}

// ----------------------------------------------------------------------
// mrc_ddc_add_value

void
mrc_ddc_amr_add_value(struct mrc_ddc *ddc,
		      int row_patch, int rowm, int row[3],
		      int col_patch, int colm, int col[3],
		      float val)
{
  struct mrc_ddc_amr *sub = mrc_ddc_amr(ddc);
  struct mrc_mat_mcsr *mcsr = &sub->mat;

  // WARNING, all elements for any given row must be added contiguously!

  assert(row_patch >= 0);
  assert(col_patch >= 0);

  int row_idx = (((rowm * sub->im[2] + row[2] - sub->ib[2]) *
		  sub->im[1] + row[1] - sub->ib[1]) *
		 sub->im[0] + row[0] - sub->ib[0]);
  int col_pidx = ((((col_patch *
  		     sub->im[3] + colm) *
  		    sub->im[2] + col[2] - sub->ib[2]) *
  		   sub->im[1] + col[1] - sub->ib[1]) *
  		  sub->im[0] + col[0] - sub->ib[0]);
  
  if (mcsr->nr_rows == 0 ||
      mcsr->rows[mcsr->nr_rows - 1].idx != row_idx ||
      mcsr->rows[mcsr->nr_rows - 1].patch != row_patch) {
    // start new row
    if (mcsr->nr_rows >= mcsr->nr_rows_alloced - 1) {
      mcsr->nr_rows_alloced *= 2;
      mcsr->rows = realloc(mcsr->rows, mcsr->nr_rows_alloced * sizeof(*mcsr->rows));
    }
    mcsr->rows[mcsr->nr_rows].patch = row_patch;
    mcsr->rows[mcsr->nr_rows].idx = row_idx;
    mcsr->rows[mcsr->nr_rows].first_entry = mcsr->nr_entries;
    mcsr->nr_rows++;
  }

  // if we already have an entry for this column in the current row, just add to it
  for (int i = mcsr->rows[mcsr->nr_rows - 1].first_entry; i < mcsr->nr_entries; i++) {
    if (mcsr->entries[i].pidx == col_pidx) {
      mcsr->entries[i].val += val;
      return;
    }
  }

  // otherwise, need to append a new entry
  if (mcsr->nr_entries >= mcsr->nr_entries_alloced) {
    mcsr->nr_entries_alloced *= 2;
    mcsr->entries = realloc(mcsr->entries, mcsr->nr_entries_alloced * sizeof(*mcsr->entries));
  }
  mcsr->entries[mcsr->nr_entries].pidx = col_pidx;
  mcsr->entries[mcsr->nr_entries].val = val;
  mcsr->nr_entries++;
}

// ----------------------------------------------------------------------
// mrc_ddc_amr_assemble

void
mrc_ddc_amr_assemble(struct mrc_ddc *ddc)
{
  struct mrc_ddc_amr *sub = mrc_ddc_amr(ddc);
  struct mrc_mat_mcsr *mcsr = &sub->mat;

  mcsr->rows[mcsr->nr_rows].first_entry = mcsr->nr_entries;
  mprintf("nr_rows %d nr_entries %d\n", mcsr->nr_rows, mcsr->nr_entries);
}

// ----------------------------------------------------------------------
// mrc_ddc_amr_fill_ghosts

static void
mrc_ddc_amr_fill_ghosts(struct mrc_ddc *ddc, struct mrc_fld *fld)
{
  struct mrc_ddc_amr *sub = mrc_ddc_amr(ddc);
  struct mrc_mat_mcsr *mcsr = &sub->mat;

  unsigned int size_of_patch = mrc_fld_ghost_dims(fld)[0] * mrc_fld_ghost_dims(fld)[1] * mrc_fld_ghost_dims(fld)[2] * mrc_fld_ghost_dims(fld)[3];

  if (ddc->size_of_type == sizeof(float)) {
    float *arr = fld->_arr;
    
    for (int row = 0; row < mcsr->nr_rows; row++) {
      int row_patch = mcsr->rows[row].patch;
      int row_idx = mcsr->rows[row].idx;
      float sum = 0.;
      for (int entry = mcsr->rows[row].first_entry;
	   entry < mcsr->rows[row + 1].first_entry; entry++) {
	int col_pidx = mcsr->entries[entry].pidx;
	float val = mcsr->entries[entry].val;
	sum += val * arr[col_pidx];
      }
      arr[row_patch * size_of_patch + row_idx] = sum;
    }
  } else {
    assert(0);
  }
}

// ----------------------------------------------------------------------
// mrc_ddc_amr_add_ghosts

static void
mrc_ddc_amr_add_ghosts(struct mrc_ddc *ddc, int mb, int me, void *ctx)
{
  // mb, me is meaningless here, make sure the caller knows
  assert(mb < 0 && me < 0);
}

// ----------------------------------------------------------------------
// mrc_ddc_amr_apply

void
mrc_ddc_amr_apply(struct mrc_ddc *ddc, struct mrc_fld *fld)
{
  mrc_ddc_amr_fill_ghosts(ddc, fld);
}

// ----------------------------------------------------------------------

#define VAR(x) (void *)offsetof(struct mrc_ddc_amr, x)
static struct param mrc_ddc_amr_descr[] = {
  { "sw"                     , VAR(sw)                      , PARAM_INT3(0, 0, 0)    },
  { "n_comp"                 , VAR(im[3])                   , PARAM_INT(0)           },
  {},
};
#undef VAR

// ======================================================================
// mrc_ddc_amr_ops

struct mrc_ddc_ops mrc_ddc_amr_ops = {
  .name                  = "amr",
  .size                  = sizeof(struct mrc_ddc_amr),
  .param_descr           = mrc_ddc_amr_descr,
  .setup                 = mrc_ddc_amr_setup,
  .destroy               = mrc_ddc_amr_destroy,
  .set_domain            = mrc_ddc_amr_set_domain,
  .get_domain            = mrc_ddc_amr_get_domain,
  //  .fill_ghosts           = mrc_ddc_amr_fill_ghosts,
  .add_ghosts            = mrc_ddc_amr_add_ghosts,
};

