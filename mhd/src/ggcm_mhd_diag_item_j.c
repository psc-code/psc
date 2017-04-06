
#include "ggcm_mhd_diag_item_private.h"

#include "ggcm_mhd_defs.h"
#include "ggcm_mhd_private.h"
#include "ggcm_mhd_diag_private.h"

#include <mrc_domain.h>

#include <stdio.h>
#include <assert.h>


// ======================================================================
// ggcm_mhd_diag_item subclass "j"

// ----------------------------------------------------------------------
// ggcm_mhd_diag_item_j_run

static void
ggcm_mhd_diag_item_j_run(struct ggcm_mhd_diag_item *item,
			 struct mrc_io *io, struct mrc_fld *f,
			 int diag_type, float plane)
{
  struct ggcm_mhd *mhd = item->diag->mhd;
  struct mrc_fld *fld_J = mrc_domain_fld_create(mhd->domain, SW_2, "jx:jy:jz");  
  mrc_fld_set_name(fld_J, "j");
  mrc_fld_setup(fld_J); 
  //currcc_f(mhd);
  ggcm_mhd_calc_currcc(mhd, mhd->fld, BX, fld_J);
  float scale_cc = mhd->ccnorm;
  ggcm_mhd_diag_c_write_one_field(io, fld_J, 0, "jx", scale_cc, diag_type, plane);
  ggcm_mhd_diag_c_write_one_field(io, fld_J, 1, "jy", scale_cc, diag_type, plane);
  ggcm_mhd_diag_c_write_one_field(io, fld_J, 2, "jz", scale_cc, diag_type, plane);
  mrc_fld_destroy(fld_J);
}

// ----------------------------------------------------------------------
// ggcm_mhd_diag_item subclass "j"

struct ggcm_mhd_diag_item_ops ggcm_mhd_diag_item_ops_j = {
  .name             = "j",
  .run              = ggcm_mhd_diag_item_j_run,
};
