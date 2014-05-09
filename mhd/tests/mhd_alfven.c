
#include "ggcm_mhd_defs.h"
#include "ggcm_mhd_private.h"
#include "ggcm_mhd_crds.h"
#include "ggcm_mhd_diag.h"
#include "ggcm_mhd_ic_private.h"

#include <mrc_fld.h>
#include <mrc_fld_as_float.h>
#include <mrc_domain.h>

#include <math.h>
#include <string.h>
#include <stdio.h>
#include <stdlib.h> 
#include <assert.h>

// ======================================================================
// ggcm_mhd_ic subclass "alfven"

struct ggcm_mhd_ic_alfven {
  float pert; // pertubation amplitude
  float B0; // initial density 
  float n0; // n0 
  float p0; // initial pressure 
  const char* pdim; 
};

// ----------------------------------------------------------------------
// ggcm_mhd_ic_alfven_run

static void
ggcm_mhd_ic_alfven_run(struct ggcm_mhd_ic *ic)
{
  struct ggcm_mhd_ic_alfven *sub = mrc_to_subobj(ic, struct ggcm_mhd_ic_alfven);
  struct ggcm_mhd *mhd = ic->mhd;
  struct mrc_fld *fld = mrc_fld_get_as(mhd->fld, FLD_TYPE);
  struct mrc_crds *crds = mrc_domain_get_crds(mhd->domain);  

  double xl[3], xh[3];
  float r[3];
  mrc_crds_get_param_double3(crds, "l", xl);
  mrc_crds_get_param_double3(crds, "h", xh);
  mrc_fld_foreach(fld, ix, iy, iz, 2, 2) {
    r[0] = MRC_CRD(crds, 0, ix);
    r[1] = MRC_CRD(crds, 1, iy);
    r[2] = MRC_CRD(crds, 2, iz);
    
    if(strcmp(sub->pdim, "x") == 0){
      float k = 2. * M_PI / (xh[0] - xl[0]);
      RR1(fld, ix,iy,iz) = sub->n0;
      PP1(fld, ix,iy,iz) = RR1(fld, ix,iy,iz);
      V1X(fld, ix,iy,iz) = 0.0;
      V1Y(fld, ix,iy,iz) = sub->pert*sin(k*r[0]);
      V1Z(fld, ix,iy,iz) = 0.0;
      B1X(fld, ix,iy,iz) = sub->B0; 
      B1Y(fld, ix,iy,iz) = sub->pert*sin(k*r[0]); 
      B1Z(fld, ix,iy,iz) = 0.0;

    } else if(strcmp(sub->pdim, "y") == 0){
      float k = 2. * M_PI / (xh[1] - xl[1]);
      RR1(fld, ix,iy,iz) = sub->n0;
      PP1(fld, ix,iy,iz) = RR1(fld, ix,iy,iz);
      V1Y(fld, ix,iy,iz) = 0.0;
      V1X(fld, ix,iy,iz) = sub->pert*sin(k*r[1]);
      V1Z(fld, ix,iy,iz) = 0.0;
      B1Y(fld, ix,iy,iz) = sub->B0; 
      B1X(fld, ix,iy,iz) = sub->pert*sin(k*r[1]); 
      B1Z(fld, ix,iy,iz) = 0.0;

    } else if(strcmp(sub->pdim, "z") == 0){
      float k = 2. * M_PI / (xh[2] - xl[2]);
      RR1(fld, ix,iy,iz) = sub->n0;
      PP1(fld, ix,iy,iz) = RR1(fld, ix,iy,iz);
      V1Y(fld, ix,iy,iz) = 0.0;
      V1X(fld, ix,iy,iz) = sub->pert*sin(k*r[2]);
      V1Z(fld, ix,iy,iz) = 0.0;
      B1Z(fld, ix,iy,iz) = sub->B0; 
      B1X(fld, ix,iy,iz) = sub->pert*sin(k*r[2]); 
      B1Y(fld, ix,iy,iz) = 0.0;
    } else {           
      assert(0); /* unknown initial condition */
    }
  } mrc_fld_foreach_end;
  
  mrc_fld_put_as(fld, mhd->fld);

  ggcm_mhd_convert_from_primitive(mhd, mhd->fld);
}

// ----------------------------------------------------------------------
// ggcm_mhd_ic_alfven_descr

#define VAR(x) (void *)offsetof(struct ggcm_mhd_ic_alfven, x)
static struct param ggcm_mhd_ic_alfven_descr[] = {
  { "pert"           , VAR(pert)          , PARAM_FLOAT(1e-3)    },
  { "B0"             , VAR(B0)            , PARAM_FLOAT(1.0)     },
  { "n0"             , VAR(n0)            , PARAM_FLOAT(1.0)     },
  { "p0"             , VAR(p0)            , PARAM_FLOAT(1.0)     }, 
  { "pdim"           , VAR(pdim)          , PARAM_STRING("x")    },  
  {},
};
#undef VAR

// ----------------------------------------------------------------------
// ggcm_mhd_ic_alfven_ops

struct ggcm_mhd_ic_ops ggcm_mhd_ic_alfven_ops = {
  .name        = "alfven",
  .size        = sizeof(struct ggcm_mhd_ic_alfven),
  .param_descr = ggcm_mhd_ic_alfven_descr,
  .run         = ggcm_mhd_ic_alfven_run,
};


// ======================================================================
// ggcm_mhd class "alfven"

// ----------------------------------------------------------------------
// ggcm_mhd_alfven_create

static void
ggcm_mhd_alfven_create(struct ggcm_mhd *mhd)
{
  ggcm_mhd_default_box(mhd);

  /* set defaults for coord arrays */
  struct mrc_crds *crds = mrc_domain_get_crds(mhd->domain);
  mrc_crds_set_type(crds, "uniform");
  mrc_crds_set_param_int(crds, "sw", SW_2);   // 'stencil width' 
  mrc_crds_set_param_double3(crds, "l", (double[3]) {  0.0, 0.0, 0.0 });
  mrc_crds_set_param_double3(crds, "h", (double[3]) {  1.0, 0.05, 0.05 });
}

static struct ggcm_mhd_ops ggcm_mhd_alfven_ops = {
  .name             = "alfven",
  .create           = ggcm_mhd_alfven_create,
};

// ======================================================================

extern struct ggcm_mhd_diag_ops ggcm_mhd_diag_c_ops;

int
main(int argc, char **argv)
{
  mrc_class_register_subclass(&mrc_class_ggcm_mhd, &ggcm_mhd_alfven_ops);  
  mrc_class_register_subclass(&mrc_class_ggcm_mhd_diag, &ggcm_mhd_diag_c_ops);
  mrc_class_register_subclass(&mrc_class_ggcm_mhd_ic, &ggcm_mhd_ic_alfven_ops);  
 
  return ggcm_mhd_main(&argc, &argv);
}

