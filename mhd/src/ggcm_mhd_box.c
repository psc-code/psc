
#include "ggcm_mhd_private.h"

#include "ggcm_mhd_crds_private.h"

#include <mrc_domain.h>

// ======================================================================
// ggcm_mhd subclass "box"
//
// this subclass of ggcm_mhd is for tests and other MHD runs, that
// use OpenGGCM/libmrc infrastructure (MHD solver, I/O), but don't
// want actual solar wind / ionosphere / ..., they rather just run
// in a (periodic by default) box.

// ----------------------------------------------------------------------
// ggcm_mhd_box_create

static void
ggcm_mhd_box_create(struct ggcm_mhd *mhd)
{
  mrc_domain_set_param_int(mhd->domain, "bcx", BC_PERIODIC);
  mrc_domain_set_param_int(mhd->domain, "bcy", BC_PERIODIC);
  mrc_domain_set_param_int(mhd->domain, "bcz", BC_PERIODIC);

  struct mrc_crds *crds = mrc_domain_get_crds(mhd->domain);
  mrc_crds_set_type(crds, "uniform");
  mrc_crds_set_param_int(crds, "sw", 2);

  ggcm_mhd_set_param_float(mhd, "isphere", 0.);
  ggcm_mhd_set_param_float(mhd, "diffsphere", 0.);
  ggcm_mhd_set_param_float(mhd, "speedlimit", 1e9);
}

// ----------------------------------------------------------------------
// ggcm_mhd_ops "box"

struct ggcm_mhd_ops ggcm_mhd_ops_box = {
  .name      = "box",
  .create    = ggcm_mhd_box_create,
};

