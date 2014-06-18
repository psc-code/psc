
#include "ggcm_mhd_ic_private.h"

#include "ggcm_mhd_defs.h"
#include "ggcm_mhd_private.h"
#include "ggcm_mhd_diag.h"

#include <mrc_domain.h>
#include <mrc_fld_as_float.h>
#include <mrc_crds.h>
#include <math.h>
#include <string.h>
#include <assert.h>

// ======================================================================
// ggcm_mhd_ic subclass "ici"

struct ggcm_mhd_ic_ici {
  float v0; // initial velocity  
  float n0; // initial density 
};
// ----------------------------------------------------------------------
// ggcm_mhd_ic_ici_run

static void
ggcm_mhd_ic_ici_run(struct ggcm_mhd_ic *ic)
{
  struct ggcm_mhd_ic_ici *sub = mrc_to_subobj(ic, struct ggcm_mhd_ic_ici);
  struct ggcm_mhd *mhd = ic->mhd;  
  struct mrc_fld *fld = mrc_fld_get_as(mhd->fld, FLD_TYPE);
  struct mrc_crds *crds = mrc_domain_get_crds(mhd->domain);  
  double xl[3], xh[3];
  float L[3], r[3];
  mrc_crds_get_param_double3(crds, "l", xl);
  mrc_crds_get_param_double3(crds, "h", xh);
  for(int i=0; i<3; i++){
    L[i] = xh[i] - xl[i];
  }

  mrc_fld_foreach(fld, ix, iy, iz, 1, 1) {
    r[0] = MRC_CRD(crds, 0, ix);
    r[1] = MRC_CRD(crds, 1, iy);
    r[2] = MRC_CRD(crds, 2, iz);
  
   // island coalescence instability 
   // based on Sullivan, Bhattacharjee & Huang 2009
    float kx = 2.0*M_PI / L[0], ky =  2.0*M_PI / L[1];
    BX(fld, ix, iy, iz) = cos(ky*r[1])*sin(kx*r[0]);
    BY(fld, ix, iy, iz) = -cos(kx*r[0])*sin(ky*r[1]); 
    // FIXME!!! I bet the 2nd should be B1Y
    RR(fld, ix, iy, iz) = sub->n0 +   
      0.5 * (1.0 - sqrt(sqr( BX(fld, ix, iy, iz))
			+ sqr( BX(fld, ix, iy, iz))));
    PP(fld, ix, iy, iz) = RR(fld, ix, iy, iz);
    VX(fld, ix, iy, iz) = sub->v0*sin(ky*r[1]);
    VY(fld, ix, iy, iz) = sub->v0*sin(kx*r[0]);
  } mrc_fld_foreach_end;

  mrc_fld_put_as(fld, mhd->fld);

  ggcm_mhd_convert_from_primitive(mhd, mhd->fld);
}

// ----------------------------------------------------------------------
// ggcm_mhd_ic_ici_descr

#define VAR(x) (void *)offsetof(struct ggcm_mhd_ic_ici, x)
static struct param ggcm_mhd_ic_ici_descr[] = {
  {"v0", VAR(v0), PARAM_FLOAT(0.1)},
  {"n0", VAR(n0), PARAM_FLOAT(1.0)},
  {},
};
#undef VAR

// ----------------------------------------------------------------------
// ggcm_mhd_ic_ici_ops

struct ggcm_mhd_ic_ops ggcm_mhd_ic_ici_ops = {
  .name        = "ici",
  .size        = sizeof(struct ggcm_mhd_ic_ici),
  .param_descr = ggcm_mhd_ic_ici_descr,
  .run         = ggcm_mhd_ic_ici_run,
};


// ======================================================================
// ggcm_mhd class "ici"

// ----------------------------------------------------------------------
// ggcm_mhd_ici_create

static void
ggcm_mhd_ici_create(struct ggcm_mhd *mhd)
{
  ggcm_mhd_default_box(mhd);
}

static struct ggcm_mhd_ops ggcm_mhd_ici_ops = {
  .name             = "ici",
  .create           = ggcm_mhd_ici_create,
};

// ======================================================================
// main

extern struct ggcm_mhd_diag_ops ggcm_mhd_diag_c_ops;

int
main(int argc, char **argv)
{
  mrc_class_register_subclass(&mrc_class_ggcm_mhd, &ggcm_mhd_ici_ops);  
  mrc_class_register_subclass(&mrc_class_ggcm_mhd_diag, &ggcm_mhd_diag_c_ops);
  mrc_class_register_subclass(&mrc_class_ggcm_mhd_ic, &ggcm_mhd_ic_ici_ops);  
 
  return ggcm_mhd_main(&argc, &argv);
}
