
#include "ggcm_mhd_defs.h"
#include "ggcm_mhd_private.h"
#include "ggcm_mhd_bnd.h"
#include "ggcm_mhd_diag.h"
#include "ggcm_mhd_ic_private.h"

#include <mrc_domain.h>
#include <mrc_crds.h>

// ======================================================================
// ggcm_mhd_ic subclass "sod"

struct ggcm_mhd_ic_sod {
  double rr_l, rr_r; // initial density 
  double pp_l, pp_r; // initial inside  pressure
};
// ----------------------------------------------------------------------
// ggcm_mhd_ic_sod_primitive

static double
ggcm_mhd_ic_sod_primitive(struct ggcm_mhd_ic *ic, int m, double crd[3])
{
  struct ggcm_mhd_ic_sod *sub = mrc_to_subobj(ic, struct ggcm_mhd_ic_sod);

  switch (m) {
  case RR: return crd[0] < .5 ? sub->rr_l : sub->rr_r;
  case PP: return crd[0] < .5 ? sub->pp_l : sub->pp_r;
  default: return 0.;
  }
}


// ----------------------------------------------------------------------
// ggcm_mhd_ic_sod_descr

#define VAR(x) (void *)offsetof(struct ggcm_mhd_ic_sod, x)
static struct param ggcm_mhd_ic_sod_descr[] = {
  { "rr_l"        , VAR(rr_l)         , PARAM_DOUBLE(1.)    },
  { "rr_r"        , VAR(rr_r)         , PARAM_DOUBLE(.125)  },
  { "pp_l"        , VAR(pp_l)         , PARAM_DOUBLE(1.0)   },
  { "pp_r"        , VAR(pp_r)         , PARAM_DOUBLE(0.1)   },
  {},
};
#undef VAR

// ----------------------------------------------------------------------
// ggcm_mhd_ic_sod_ops

struct ggcm_mhd_ic_ops ggcm_mhd_ic_sod_ops = {
  .name        = "sod",
  .size        = sizeof(struct ggcm_mhd_ic_sod),
  .param_descr = ggcm_mhd_ic_sod_descr,
  .primitive   = ggcm_mhd_ic_sod_primitive,
};

// ======================================================================
// ggcm_mhd class "sod"

// ----------------------------------------------------------------------
// ggcm_mhd_sod_create

static void
ggcm_mhd_sod_create(struct ggcm_mhd *mhd)
{
  ggcm_mhd_default_box(mhd);

  ggcm_mhd_bnd_set_type(mhd->bnd, "open_x");
  mrc_domain_set_param_int(mhd->domain, "bcx", BC_NONE);

  mrc_domain_set_param_int3(mhd->domain, "m", (int [3]) { 128, 1, 1 });

  // set defaults for coord arrays
  struct mrc_crds *crds = mrc_domain_get_crds(mhd->domain);
  mrc_crds_set_type(crds, "uniform");
  mrc_crds_set_param_double3(crds, "l", (double[3]) {  0.0, 0.0, 0.0 });
  mrc_crds_set_param_double3(crds, "h", (double[3]) {  1.0, 0.1, 0.1 });
}

static struct ggcm_mhd_ops ggcm_mhd_sod_ops = {
  .name             = "sod",
  .create           = ggcm_mhd_sod_create,
};

// ======================================================================

extern struct ggcm_mhd_diag_ops ggcm_mhd_diag_c_ops;

int
main(int argc, char **argv)
{
  mrc_class_register_subclass(&mrc_class_ggcm_mhd, &ggcm_mhd_sod_ops);  
  mrc_class_register_subclass(&mrc_class_ggcm_mhd_diag, &ggcm_mhd_diag_c_ops);
  mrc_class_register_subclass(&mrc_class_ggcm_mhd_ic, &ggcm_mhd_ic_sod_ops);  
 
  return ggcm_mhd_main(&argc, &argv);
}

