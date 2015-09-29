
#include "ggcm_mhd_diag_item_private.h"

#include "ggcm_mhd_defs.h"
#include "ggcm_mhd_private.h"
#include "ggcm_mhd_diag_private.h"

#include <mrc_domain.h>
#include <mrc_fld_as_float.h>

#include <stdio.h>
#include <assert.h>

// ======================================================================
// ggcm_mhd_diag_item subclass "rank"

// ----------------------------------------------------------------------
// ggcm_mhd_diag_item_rank_run

static void
ggcm_mhd_diag_item_rank_run(struct ggcm_mhd_diag_item *item,
                            struct mrc_io *io, struct mrc_fld *fld,
                            int diag_type, float plane)
{
  int bnd = fld->_nr_ghosts;

  struct ggcm_mhd *mhd = item->diag->mhd;
  struct mrc_fld *fld_rank = mrc_domain_fld_create(mhd->domain, bnd, "rank");
  mrc_fld_set_type(fld_rank, FLD_TYPE);
  mrc_fld_setup(fld_rank);
  
  // ggcm_mhd_calc_rank(mhd, fld, fld_rank);
  for (int p = 0; p < mrc_fld_nr_patches(fld); p++) {
    struct mrc_patch_info info;
    mrc_domain_get_local_patch_info(fld->_domain, p, &info);
    mrc_fld_foreach(fld, ix,iy,iz, bnd, bnd) {
      M3(fld_rank, 0, ix,iy,iz, p) = info.rank;
    } mrc_fld_foreach_end;
  }
  ggcm_mhd_diag_c_write_one_fld(io, fld_rank, diag_type, plane);

  mrc_fld_destroy(fld_rank);
}

// ----------------------------------------------------------------------
// ggcm_mhd_diag_item subclass "rank"

struct ggcm_mhd_diag_item_ops ggcm_mhd_diag_item_ops_rank = {
  .name             = "rank",
  .run              = ggcm_mhd_diag_item_rank_run,
};
