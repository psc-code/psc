
#include "ggcm_mhd_bnd_private.h"

// ======================================================================
// ggcm_mhd_bnd subclass "none"

// ----------------------------------------------------------------------
// ggcm_mhd_bnd_none_fill_ghosts

static void
ggcm_mhd_bnd_none_fill_ghosts(struct ggcm_mhd_bnd *bnd, struct mrc_fld *fld,
			      float bntim)
{
}

// ----------------------------------------------------------------------
// ggcm_mhd_bnd subclass "none"

struct ggcm_mhd_bnd_ops ggcm_mhd_bnd_ops_none = {
  .name             = "none",
  .fill_ghosts      = ggcm_mhd_bnd_none_fill_ghosts,
};

