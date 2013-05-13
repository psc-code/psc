
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
			    struct mrc_io *io, struct mrc_f3 *f,
			    int diag_type, float plane)
{
  struct ggcm_mhd *mhd = item->diag->mhd;
  struct mrc_f3 *fld_divB = mrc_domain_f3_create(mhd->domain, SW_2);
  mrc_f3_setup(fld_divB);
  mrc_f3_set_comp_name(fld_divB, 0, "divB");
  
  ggcm_mhd_calc_divb(mhd, mhd->flds_base, fld_divB);
  ggcm_mhd_diag_c_write_one_f3(io, fld_divB, diag_type, plane);

  mrc_f3_destroy(fld_divB);

}

// ----------------------------------------------------------------------
// ggcm_mhd_diag_item subclass "divb"

struct ggcm_mhd_diag_item_ops ggcm_mhd_diag_item_ops_divb = {
  .name             = "divb",
  .run              = ggcm_mhd_diag_item_divb_run,
};

