
#include "ggcm_mhd_defs.h"
#include "ggcm_mhd_private.h"
#include "ggcm_mhd_bnd.h"
#include "ggcm_mhd_diag.h"
#include "ggcm_mhd_ic_private.h"

#include <mrc_domain.h>
#include <mrc_crds.h>

// ======================================================================
// ggcm_mhd_ic subclass "shocktube"

struct ggcm_mhd_ic_shocktube {
  double rr_l, rr_r; // initial density 
  double pp_l, pp_r; // initial inside  pressure
  double vx_l, vx_r; // initial vx
  double vy_l, vy_r; // initial vy
  double vz_l, vz_r; // initial vz
  double bx;         // initial Bx
  double by_l, by_r; // initial By
  double bz_l, bz_r; // initial Bz
};
// ----------------------------------------------------------------------
// ggcm_mhd_ic_shocktube_primitive

static double
ggcm_mhd_ic_shocktube_primitive(struct ggcm_mhd_ic *ic, int m, double crd[3])
{
  struct ggcm_mhd_ic_shocktube *sub = mrc_to_subobj(ic, struct ggcm_mhd_ic_shocktube);

  switch (m) {
  case RR: return crd[0] < .5 ? sub->rr_l : sub->rr_r;
  case PP: return crd[0] < .5 ? sub->pp_l : sub->pp_r;
  case VX: return crd[0] < .5 ? sub->vx_l : sub->vx_r;
  case VY: return crd[0] < .5 ? sub->vy_l : sub->vy_r;
  case VZ: return crd[0] < .5 ? sub->vz_l : sub->vz_r;
  case BX: return sub->bx;
  case BY: return crd[0] < .5 ? sub->by_l : sub->by_r;
  case BZ: return crd[0] < .5 ? sub->bz_l : sub->bz_r;
  default: return 0.;
  }
}


// ----------------------------------------------------------------------
// ggcm_mhd_ic_shocktube_descr

#define VAR(x) (void *)offsetof(struct ggcm_mhd_ic_shocktube, x)
static struct param ggcm_mhd_ic_shocktube_descr[] = {
  { "rr_l"        , VAR(rr_l)         , PARAM_DOUBLE(1.)    },
  { "rr_r"        , VAR(rr_r)         , PARAM_DOUBLE(.125)  },
  { "pp_l"        , VAR(pp_l)         , PARAM_DOUBLE(1.)    },
  { "pp_r"        , VAR(pp_r)         , PARAM_DOUBLE(.1)    },
  { "vx_l"        , VAR(vx_l)         , PARAM_DOUBLE(0.)    },
  { "vx_r"        , VAR(vx_r)         , PARAM_DOUBLE(0.)    },
  { "vy_l"        , VAR(vy_l)         , PARAM_DOUBLE(0.)    },
  { "vy_r"        , VAR(vy_r)         , PARAM_DOUBLE(0.)    },
  { "vz_l"        , VAR(vz_l)         , PARAM_DOUBLE(0.)    },
  { "vz_r"        , VAR(vz_r)         , PARAM_DOUBLE(0.)    },
  { "bx"          , VAR(bx)           , PARAM_DOUBLE(.75)   },
  { "by_l"        , VAR(by_l)         , PARAM_DOUBLE(1.)    },
  { "by_r"        , VAR(by_r)         , PARAM_DOUBLE(-1.)   },
  { "bz_l"        , VAR(bz_l)         , PARAM_DOUBLE(0.)    },
  { "bz_r"        , VAR(bz_r)         , PARAM_DOUBLE(0.)    },
  {},
};
#undef VAR

// ----------------------------------------------------------------------
// ggcm_mhd_ic_shocktube_ops

struct ggcm_mhd_ic_ops ggcm_mhd_ic_shocktube_ops = {
  .name        = "shocktube",
  .size        = sizeof(struct ggcm_mhd_ic_shocktube),
  .param_descr = ggcm_mhd_ic_shocktube_descr,
  .primitive   = ggcm_mhd_ic_shocktube_primitive,
};

// ======================================================================
// ggcm_mhd class "shocktube"

// ----------------------------------------------------------------------
// ggcm_mhd_shocktube_create

static void
ggcm_mhd_shocktube_create(struct ggcm_mhd *mhd)
{
  ggcm_mhd_default_box(mhd);
  mhd->par.gamm = 2.;

  ggcm_mhd_bnd_set_type(mhd->bnd, "open_x");
  mrc_domain_set_param_int(mhd->domain, "bcx", BC_NONE);

  mrc_domain_set_param_int3(mhd->domain, "m", (int [3]) { 128, 1, 1 });

  // set defaults for coord arrays
  struct mrc_crds *crds = mrc_domain_get_crds(mhd->domain);
  mrc_crds_set_type(crds, "uniform");
  mrc_crds_set_param_double3(crds, "l", (double[3]) {  0.0, 0.0, 0.0 });
  mrc_crds_set_param_double3(crds, "h", (double[3]) {  1.0, 0.1, 0.1 });
}

static struct ggcm_mhd_ops ggcm_mhd_shocktube_ops = {
  .name             = "shocktube",
  .create           = ggcm_mhd_shocktube_create,
};

// ======================================================================

extern struct ggcm_mhd_diag_ops ggcm_mhd_diag_c_ops;

int
main(int argc, char **argv)
{
  mrc_class_register_subclass(&mrc_class_ggcm_mhd, &ggcm_mhd_shocktube_ops);  
  mrc_class_register_subclass(&mrc_class_ggcm_mhd_diag, &ggcm_mhd_diag_c_ops);
  mrc_class_register_subclass(&mrc_class_ggcm_mhd_ic, &ggcm_mhd_ic_shocktube_ops);  
 
  return ggcm_mhd_main(&argc, &argv);
}

