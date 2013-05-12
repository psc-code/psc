
#include "ggcm_mhd_crds_private.h"

#include "ggcm_mhd_defs.h"
#include "ggcm_mhd_private.h"

#include <assert.h>

// ======================================================================
// ggcm_mhd_crds type "c"

// ----------------------------------------------------------------------
// ggcm_mhd_crds_c_setup

static void
ggcm_mhd_crds_c_setup(struct ggcm_mhd_crds *crds)
{
  assert(crds->domain);
  struct mrc_patch_info info;
  mrc_domain_get_local_patch_info(crds->domain, 0, &info);

  for (int d = 0; d < 3; d++) {
    mrc_f1_set_param_int(crds->f1[d], "offx", -BND);
    mrc_f1_set_param_int(crds->f1[d], "dimsx", info.ldims[d] + 2*BND);
  }

  ggcm_mhd_crds_setup_super(crds);
}

// ----------------------------------------------------------------------
// ggcm_mhd_crds subclass "c"

struct ggcm_mhd_crds_ops ggcm_mhd_crds_c_ops = {
  .name       = "c",
  .setup      = ggcm_mhd_crds_c_setup,
};
