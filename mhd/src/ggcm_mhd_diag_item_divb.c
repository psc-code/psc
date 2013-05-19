
#include "ggcm_mhd_diag_item_private.h"

#include "ggcm_mhd_defs.h"
#include "ggcm_mhd_private.h"
#include "ggcm_mhd_diag_private.h"

#include <mrc_domain.h>

#include <stdio.h>
#include <assert.h>

// ======================================================================
// ggcm_mhd_diag_item subclass "divb"

// ----------------------------------------------------------------------
// ggcm_mhd_diag_item_divb_run

static void
ggcm_mhd_diag_item_divb_run(struct ggcm_mhd_diag_item *item,
			    struct mrc_io *io, struct mrc_fld *f,
			    int diag_type, float plane)
{
  struct ggcm_mhd *mhd = item->diag->mhd;
  struct mrc_fld *fld_divB = mrc_domain_fld_create(mhd->domain, SW_2, "divB");
  mrc_fld_setup(fld_divB);
  
  ggcm_mhd_calc_divb(mhd, mhd->flds_base, fld_divB);
  ggcm_mhd_diag_c_write_one_fld(io, fld_divB, diag_type, plane);

  mrc_fld_destroy(fld_divB);
}

// ----------------------------------------------------------------------
// ggcm_mhd_diag_item subclass "divb"

struct ggcm_mhd_diag_item_ops ggcm_mhd_diag_item_ops_divb = {
  .name             = "divb",
  .run              = ggcm_mhd_diag_item_divb_run,
};

