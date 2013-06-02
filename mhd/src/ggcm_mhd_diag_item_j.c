
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
			 struct mrc_io *io, struct mrc_f3 *f,
			 int diag_type, float plane)
{
  struct ggcm_mhd *mhd = item->diag->mhd;
  struct mrc_f3 *fld_J = mrc_domain_f3_create(mhd->domain, 3);  
  mrc_f3_set_name(fld_J, "j");
  mrc_f3_set_nr_comps(fld_J, 3); 
  mrc_f3_set_comp_name(fld_J, 0, "jx");
  mrc_f3_set_comp_name(fld_J, 1, "jy");
  mrc_f3_set_comp_name(fld_J, 2, "jz");
  mrc_f3_setup(fld_J); 
  //currcc_f(mhd);
  ggcm_mhd_calc_currcc(mhd, mhd->flds_base, _B1X, fld_J);
  //ggcm_mhd_calc_curl(mhd, mhd->flds_base, _B1X, fld_J, tmp_sc);
  // FIXME float scale_cc = mhd->par.ccnorm;
  ggcm_mhd_diag_c_write_one_f3(io, fld_J, diag_type, plane);  
  mrc_f3_destroy(fld_J);
}

// ----------------------------------------------------------------------
// ggcm_mhd_diag_item subclass "j"

struct ggcm_mhd_diag_item_ops ggcm_mhd_diag_item_ops_j = {
  .name             = "j",
  .run              = ggcm_mhd_diag_item_j_run,
};
