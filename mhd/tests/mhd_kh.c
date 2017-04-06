
#include "ggcm_mhd_defs.h"
#include "ggcm_mhd_private.h"
#include "ggcm_mhd_crds.h"
#include "ggcm_mhd_diag.h"
#include "ggcm_mhd_ic_private.h"

#include <mrc_fld.h>
#include <mrc_fld_as_double.h>
#include <mrc_domain.h>
#include <mrc_crds.h> 

#include <math.h>
#include <string.h>
#include <stdio.h>
#include <stdlib.h> 
#include <assert.h>

// ----------------------------------------------------------------------

static double
random_double()
{
  return (double) random() / RAND_MAX;
}

// ======================================================================
// ggcm_mhd_ic subclass "kh"

struct ggcm_mhd_ic_kh {
  double rho0; // initial density 0 
  double rho1; // initial density 1  
  double v0; // velocity 0 
  double v1; // velocity 1
  double p0; // initial pressure
  double B0x; // initial Bx
  double B0y; // initial By
  double B0z; // initial Bz
  double B0z_harris; // initial Bz harris
  double pert; // initial pertubation amplitude
  double pert_random; // random pertubation amplitude
  double lambda; // shear layer width
  double sigma; // width of perturbation
};

// ----------------------------------------------------------------------
// ggcm_mhd_ic_kh_descr

#define VAR(x) (void *)offsetof(struct ggcm_mhd_ic_kh, x)
static struct param ggcm_mhd_ic_kh_descr[] = {
  { "rho0"          , VAR(rho0)          , PARAM_DOUBLE(2.)      },
  { "rho1"          , VAR(rho1)          , PARAM_DOUBLE(1.)      },
  { "v0"            , VAR(v0)            , PARAM_DOUBLE(.5)      },
  { "v1"            , VAR(v1)            , PARAM_DOUBLE(-.5)     },
  { "B0x"           , VAR(B0x)           , PARAM_DOUBLE(0.)      },
  { "B0y"           , VAR(B0y)           , PARAM_DOUBLE(0.)      },
  { "B0z"           , VAR(B0z)           , PARAM_DOUBLE(0.)      },
  { "B0z_harris"    , VAR(B0z_harris)    , PARAM_DOUBLE(0.)      },
  { "p0"            , VAR(p0)            , PARAM_DOUBLE(2.5)     },
  { "pert"          , VAR(pert)          , PARAM_DOUBLE(1e-2)    },
  { "pert_random"   , VAR(pert_random)   , PARAM_DOUBLE(0.)      },
  { "lambda"        , VAR(lambda)        , PARAM_DOUBLE(.05)     },
  { "sigma"         , VAR(sigma)         , PARAM_DOUBLE(.2)      },
  {},
};
#undef VAR

// ----------------------------------------------------------------------
// ggcm_mhd_ic_kh0_primitive
//
// single mode and random perturbation, finite width shear

static double
ggcm_mhd_ic_kh0_primitive(struct ggcm_mhd_ic *ic, int m, double crd[3])
{
  struct ggcm_mhd_ic_kh *sub = mrc_to_subobj(ic, struct ggcm_mhd_ic_kh);
  double xx = crd[0], yy = crd[1];

  double s = 1. + .5 * (tanh((yy - .25) / sub->lambda) - tanh((yy + .25) / sub->lambda));
  
  switch (m) {
  case RR: return sub->rho0 * s + sub->rho1 * (1. - s);

  // VX: shear flow + random perturbation
  case VX: return sub->v0 * s + sub->v1 * (1. - s)
      + sub->pert_random * (random_double() - .5);
    
  // BY: single mode + random perturbation
  case VY: return sub->pert * sin(2.*M_PI * xx) * exp(-sqr(yy) / sqr(sub->sigma))
      + sub->pert_random * (random_double() - .5);
  
  case PP: return sub->p0 
      + .5 * sqr(sub->B0z_harris) / sqr(cosh((yy - .25) / sub->lambda))
      + .5 * sqr(sub->B0z_harris) / sqr(cosh((yy + .25) / sub->lambda));

  case BX: return sub->B0x; 
  case BY: return sub->B0y;
  case BZ: return sub->B0z + sub->B0z_harris * s - sub->B0z_harris * (1. - s);
  default: return 0.;
  }
}

// ----------------------------------------------------------------------
// ggcm_mhd_ic_kh0_ops

struct ggcm_mhd_ic_ops ggcm_mhd_ic_kh0_ops = {
  .name        = "kh",
  .size        = sizeof(struct ggcm_mhd_ic_kh),
  .param_descr = ggcm_mhd_ic_kh_descr,
  .primitive   = ggcm_mhd_ic_kh0_primitive,
};



// ----------------------------------------------------------------------
// ggcm_mhd_ic_kh1_primitive
//
// athena-like: random perturbations

static double
ggcm_mhd_ic_kh1_primitive(struct ggcm_mhd_ic *ic, int m, double crd[3])
{
  struct ggcm_mhd_ic_kh *sub = mrc_to_subobj(ic, struct ggcm_mhd_ic_kh);
  double yy = crd[1];

  switch (m) {
  case RR: return (fabs(yy) < .25 ? sub->rho0 : sub->rho1);
  case VX: return (fabs(yy) < .25 ? sub->v1   : sub->v0) 
      + sub->pert_random * (random_double() - .5);
  case VY: return 0.f
      + sub->pert_random * (random_double() - .5);
  case PP: return sub->p0;
  case BX: return sub->B0x;
  case BY: return sub->B0y;
  case BZ: return sub->B0z;
  default: return 0.f;
  }
}

// ----------------------------------------------------------------------
// ggcm_mhd_ic_kh1_ops

struct ggcm_mhd_ic_ops ggcm_mhd_ic_kh1_ops = {
  .name        = "kh",
  .size        = sizeof(struct ggcm_mhd_ic_kh),
  .param_descr = ggcm_mhd_ic_kh_descr,
  .primitive   = ggcm_mhd_ic_kh1_primitive,
};



// ======================================================================
// ggcm_mhd class "kh"

// ----------------------------------------------------------------------
// ggcm_mhd_kh_create

static void
ggcm_mhd_kh_create(struct ggcm_mhd *mhd)
{
  ggcm_mhd_default_box(mhd);

  // default mesh size
  mrc_domain_set_param_int3(mhd->domain, "m", (int[3]) { 64, 64, 1 });
  
  // default domain size
  struct mrc_crds *crds = mrc_domain_get_crds(mhd->domain);
  mrc_crds_set_type(crds, "uniform");
  mrc_crds_set_param_int(crds, "sw", SW_2);   // 'stencil width' 
  mrc_crds_set_param_double3(crds, "l", (double[3]) {  0.0, 0.0, 0.0 });
  mrc_crds_set_param_double3(crds, "h", (double[3]) {  1.0, 1.0, 0.1 });
}

// ----------------------------------------------------------------------
// ggcm_mhd_kh_ops 

static struct ggcm_mhd_ops ggcm_mhd_kh_ops = {
  .name             = "kh",
  .create           = ggcm_mhd_kh_create,
};

// ======================================================================

extern struct ggcm_mhd_diag_ops ggcm_mhd_diag_c_ops;

int
main(int argc, char **argv)
{
  mrc_class_register_subclass(&mrc_class_ggcm_mhd, &ggcm_mhd_kh_ops);  
  mrc_class_register_subclass(&mrc_class_ggcm_mhd_diag, &ggcm_mhd_diag_c_ops); 
  mrc_class_register_subclass(&mrc_class_ggcm_mhd_ic, &ggcm_mhd_ic_kh0_ops);  
  mrc_class_register_subclass(&mrc_class_ggcm_mhd_ic, &ggcm_mhd_ic_kh1_ops);  

  return ggcm_mhd_main(&argc, &argv);
}

