
#include "mrc_ddc_private.h"

#include <mrc_mat.h>
#include <mrc_domain.h>
#include <mrc_params.h>
#include <assert.h>
#include <stdlib.h>
#include <stdio.h>

// ======================================================================
// mrc_ddc_amr

struct mrc_ddc_amr {
  struct mrc_mat *mat;

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
  assert(sub->domain);

  int ldims[3];
  mrc_domain_get_param_int3(sub->domain, "m", ldims);
  int size = 1;
  // needs to be compatible with how mrc_fld indexes its fields
  for (int d = 0; d < 3; d++) {
    sub->ib[d] = -sub->sw[d];
    sub->im[d] = ldims[d] + 2 * sub->sw[d];
    size *= sub->im[d];
  }

  int nr_patches;
  mrc_domain_get_param_int(sub->domain, "nr_patches", &nr_patches);
  size *= sub->im[3]; // # components
  size *= nr_patches;

  sub->mat = mrc_mat_create(mrc_ddc_comm(ddc));
  mrc_mat_set_type(sub->mat, "csr_mpi");
  mprintf("size = %d %d im %d\n", size, nr_patches, sub->im[3]);
  mrc_mat_set_param_int(sub->mat, "m", size);
  mrc_mat_set_param_int(sub->mat, "n", size);
  mrc_mat_set_from_options(sub->mat); // to allow changing matrix type
  mrc_mat_setup(sub->mat);
}

// ----------------------------------------------------------------------
// mrc_ddc_amr_destroy

static void
mrc_ddc_amr_destroy(struct mrc_ddc *ddc)
{
  struct mrc_ddc_amr *sub = mrc_ddc_amr(ddc);
  mrc_mat_destroy(sub->mat);
}

// ----------------------------------------------------------------------
// mrc_ddc_add_value

void
mrc_ddc_amr_add_value(struct mrc_ddc *ddc,
		      int row_patch, int rowm, int row[3],
		      int col_patch, int colm, int col[3],
		      double val)
{
  struct mrc_ddc_amr *sub = mrc_ddc_amr(ddc);

  assert(row_patch >= 0);
  assert(row[0] >= sub->ib[0] && row[0] < sub->ib[0] + sub->im[0]);
  assert(row[1] >= sub->ib[1] && row[1] < sub->ib[1] + sub->im[1]);
  assert(row[2] >= sub->ib[2] && row[2] < sub->ib[2] + sub->im[2]);

  assert(col_patch >= 0);
  assert(col[0] >= sub->ib[0] && col[0] < sub->ib[0] + sub->im[0]);
  assert(col[1] >= sub->ib[1] && col[1] < sub->ib[1] + sub->im[1]);
  assert(col[2] >= sub->ib[2] && col[2] < sub->ib[2] + sub->im[2]);

  int row_idx = ((((row_patch *
		    sub->im[3] + rowm) *
		   sub->im[2] + row[2] - sub->ib[2]) *
		  sub->im[1] + row[1] - sub->ib[1]) *
		 sub->im[0] + row[0] - sub->ib[0]);
  int col_idx = ((((col_patch *
		    sub->im[3] + colm) *
		   sub->im[2] + col[2] - sub->ib[2]) *
		  sub->im[1] + col[1] - sub->ib[1]) *
		 sub->im[0] + col[0] - sub->ib[0]);

  mrc_mat_add_value(sub->mat, row_idx, col_idx, val);
}

// ----------------------------------------------------------------------
// mrc_ddc_amr_assemble

void
mrc_ddc_amr_assemble(struct mrc_ddc *ddc)
{
  struct mrc_ddc_amr *sub = mrc_ddc_amr(ddc);

  mrc_mat_assemble(sub->mat);
}

// ----------------------------------------------------------------------
// mrc_ddc_amr_apply

void
mrc_ddc_amr_apply(struct mrc_ddc *ddc, struct mrc_fld *fld)
{
  struct mrc_ddc_amr *sub = mrc_ddc_amr(ddc);

  mrc_mat_apply_in_place(sub->mat, fld->_nd->vec);
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
};

