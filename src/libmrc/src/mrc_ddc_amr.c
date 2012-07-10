
#include "mrc_ddc_private.h"

#include <mrc_domain.h>
#include <mrc_params.h>
#include <assert.h>
#include <stdlib.h>
#include <stdio.h>

#define mrc_ddc_amr(ddc) mrc_to_subobj(ddc, struct mrc_ddc_amr)

// ----------------------------------------------------------------------
// mrc_ddc_amr_set_domain

static void
mrc_ddc_amr_set_domain(struct mrc_ddc *ddc, struct mrc_domain *domain)
{
  struct mrc_ddc_amr *amr = mrc_ddc_amr(ddc);
  amr->domain = domain;
}

// ----------------------------------------------------------------------
// mrc_ddc_amr_get_domain

static struct mrc_domain *
mrc_ddc_amr_get_domain(struct mrc_ddc *ddc)
{
  struct mrc_ddc_amr *amr = mrc_ddc_amr(ddc);
  return amr->domain;
}

// ----------------------------------------------------------------------
// mrc_ddc_amr_setup

static void
mrc_ddc_amr_setup(struct mrc_ddc *ddc)
{
  struct mrc_ddc_amr *amr = mrc_ddc_amr(ddc);
  assert(amr->domain);

  int ldims[3];
  mrc_domain_get_param_int3(amr->domain, "m", ldims);
  // needs to be compatible with how mrc_m3 indexes its fields
  for (int d = 0; d < 3; d++) {
    amr->ib[d] = -amr->sw;
    amr->im[d] = ldims[d] + 2 * amr->sw;
  }

  amr->nr_rows_alloced = 1000;
  amr->nr_entries_alloced = 2000;

  amr->rows = calloc(amr->nr_rows_alloced, sizeof(*amr->rows));
  amr->entries = calloc(amr->nr_entries_alloced, sizeof(*amr->entries));

  amr->nr_entries = 0;
  amr->nr_rows = 0;
}

// ----------------------------------------------------------------------
// mrc_ddc_amr_destroy

static void
mrc_ddc_amr_destroy(struct mrc_ddc *ddc)
{
  struct mrc_ddc_amr *amr = mrc_ddc_amr(ddc);

  free(amr->rows);
  free(amr->entries);
}

// ----------------------------------------------------------------------
// mrc_ddc_add_value

void
mrc_ddc_amr_add_value(struct mrc_ddc *ddc,
		      int row_patch, int rowm, int row[3],
		      int col_patch, int colm, int col[3],
		      float val)
{
  struct mrc_ddc_amr *amr = mrc_ddc_amr(ddc);

  // WARNING, all elements for any given row must be added contiguously!

  assert(row_patch >= 0);
  assert(col_patch >= 0);

  int row_idx = (((rowm * amr->im[2] + row[2] - amr->ib[2]) *
		  amr->im[1] + row[1] - amr->ib[1]) *
		 amr->im[0] + row[0] - amr->ib[0]);
  int col_idx = (((colm * amr->im[2] + col[2] - amr->ib[2]) *
		  amr->im[1] + col[1] - amr->ib[1]) *
		 amr->im[0] + col[0] - amr->ib[0]);
  
  if (amr->nr_rows == 0 ||
      amr->rows[amr->nr_rows - 1].idx != row_idx ||
      amr->rows[amr->nr_rows - 1].patch != row_patch) {
    // start new row
    if (amr->nr_rows >= amr->nr_rows_alloced - 1) {
      amr->nr_rows_alloced *= 2;
      amr->rows = realloc(amr->rows, amr->nr_rows_alloced * sizeof(*amr->rows));
    }
    amr->rows[amr->nr_rows].patch = row_patch;
    amr->rows[amr->nr_rows].idx = row_idx;
    amr->rows[amr->nr_rows].first_entry = amr->nr_entries;
    amr->nr_rows++;
  }

  // if we already have an entry for this column in the current row, just add to it
  for (int i = amr->rows[amr->nr_rows - 1].first_entry; i < amr->nr_entries; i++) {
    if (amr->entries[i].patch == col_patch && amr->entries[i].idx == col_idx) {
      amr->entries[i].val += val;
      return;
    }
  }

  // otherwise, need to append a new entry
  if (amr->nr_entries >= amr->nr_entries_alloced) {
    amr->nr_entries_alloced *= 2;
    amr->entries = realloc(amr->entries, amr->nr_entries_alloced * sizeof(*amr->entries));
  }
  amr->entries[amr->nr_entries].patch = col_patch;
  amr->entries[amr->nr_entries].idx = col_idx;
  amr->entries[amr->nr_entries].val = val;
  amr->nr_entries++;
}

// ----------------------------------------------------------------------
// mrc_ddc_amr_assemble

void
mrc_ddc_amr_assemble(struct mrc_ddc *ddc)
{
  struct mrc_ddc_amr *amr = mrc_ddc_amr(ddc);

  amr->rows[amr->nr_rows].first_entry = amr->nr_entries;
  mprintf("nr_rows %d nr_entries %d\n", amr->nr_rows, amr->nr_entries);
}

// ----------------------------------------------------------------------
// mrc_ddc_amr_fill_ghosts

static void
mrc_ddc_amr_fill_ghosts(struct mrc_ddc *ddc, int mb, int me, void *ctx)
{
  float **fldp = ctx;
  struct mrc_ddc_amr *amr = mrc_ddc_amr(ddc);

  for (int row = 0; row < amr->nr_rows; row++) {
    int row_patch = amr->rows[row].patch;
    int row_idx = amr->rows[row].idx;
    float sum = 0.;
    for (int entry = amr->rows[row].first_entry;
	 entry < amr->rows[row + 1].first_entry; entry++) {
      int col_patch =  amr->entries[entry].patch;
      int col_idx = amr->entries[entry].idx;
      float val = amr->entries[entry].val;
      sum += val * fldp[col_patch][col_idx];
    }
    fldp[row_patch][row_idx] = sum;
  }
}

// ----------------------------------------------------------------------
// mrc_ddc_amr_apply

void
mrc_ddc_amr_apply(struct mrc_ddc *ddc, struct mrc_m3 *fld)
{
  float **fldp = malloc(fld->nr_patches * sizeof(*fldp));
  for (int p = 0; p < fld->nr_patches; p++) {
    fldp[p] = mrc_m3_patch_get(fld, p)->arr;
  }
  mrc_ddc_amr_fill_ghosts(ddc, -1, -1, fldp);

  free(fldp);
}

// ----------------------------------------------------------------------

#define VAR(x) (void *)offsetof(struct mrc_ddc_amr, x)
static struct param mrc_ddc_amr_descr[] = {
  { "sw"                     , VAR(sw)                      , PARAM_INT(0)           },
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
  .fill_ghosts           = mrc_ddc_amr_fill_ghosts,
};

