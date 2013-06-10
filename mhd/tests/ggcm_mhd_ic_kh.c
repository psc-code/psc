
#include "ggcm_mhd_ic_private.h"

#include "ggcm_mhd_defs.h"
#include "ggcm_mhd_private.h"
#include "ggcm_mhd_flds.h"

#include <mrc_domain.h>
#include <mrc_crds.h>
#include <math.h>
#include <string.h>
#include <assert.h>

#define ggcm_mhd_cweno(obj) mrc_to_subobj(obj, struct mhd)

// ======================================================================
// ggcm_mhd_ic subclass "kh"

struct ggcm_mhd_ic_kh {
  float mpermi;  
  float pert; // initial pertubation amplitude
  float r0; // initial density 0 
  float r1; // initial density 1  
  float v0; // velocity 0 
  float v1; // velocity 1
  float B0; // initial B
  float p0; // initial pressure
  float lambda; // wave number 
};
// ----------------------------------------------------------------------
// ggcm_mhd_ic_kh_run
// Kelvin-Helmholtz Instability test 
// default settings taken from 
// http://www.astro.princeton.edu/~jstone/Athena/tests/kh/kh.html

static void
ggcm_mhd_ic_kh_run(struct ggcm_mhd_ic *ic)
{
  struct ggcm_mhd_ic_kh *sub = mrc_to_subobj(ic, struct ggcm_mhd_ic_kh);
  struct ggcm_mhd *gmhd = ic->mhd;  
  struct mrc_fld *fld = gmhd->fld;
  struct mrc_crds *crds = mrc_domain_get_crds(gmhd->domain);  
  float xl[3], xh[3],  xmid[3], L[3], r[3];
  mrc_domain_set_param_int(gmhd->domain, "bcx", BC_PERIODIC);
  mrc_domain_set_param_int(gmhd->domain, "bcy", BC_PERIODIC);
  mrc_domain_set_param_int(gmhd->domain, "bcz", BC_PERIODIC);

  mrc_crds_get_xl_xh(crds, xl, xh);
  for(int i=0; i<3; i++){
    L[i] = xh[i] - xl[i];
    xmid[i] = 0.5 * (xh[i] + xl[i]);
  }
  float gamma = gmhd->par.gamm;
  mrc_fld_foreach(fld, ix, iy, iz, 1, 1) {
    r[0] = MRC_CRD(crds, 0, ix);
    r[1] = MRC_CRD(crds, 1, iy);
    r[2] = MRC_CRD(crds, 2, iz);

    float wave1 = sin(2.* sub->lambda * M_PI * (r[0] - xmid[0]) / L[0]); 
 
    if(fabs(r[1]) < 0.25*L[1]){
      MRC_F3(fld, _RR1, ix, iy, iz) = sub->r0;
      MRC_F3(fld, _RV1X , ix, iy, iz) = sub->r0 * sub->v0;
      MRC_F3(fld, _B1X, ix, iy, iz) = sub->B0;
      MRC_F3(fld, _RV1Y , ix, iy, iz) = sub->r0 * sub->pert * wave1;
    }else{
      MRC_F3(fld, _RR1, ix, iy, iz) = sub->r1;
      MRC_F3(fld, _RV1X , ix, iy, iz) = sub->r1 * sub->v1;
      MRC_F3(fld, _B1X, ix, iy, iz) = sub->B0;
      MRC_F3(fld, _RV1Y , ix, iy, iz) = sub->r1 * sub->pert * wave1;
    }   
    MRC_F3(fld, _RV1Z , ix, iy, iz) = 0.0;
    MRC_F3(fld, _UU1 , ix, iy, iz) = sub->p0/ 
      (gamma - 1.f) +
      .5f * (sqr(MRC_F3(fld, _RV1X, ix, iy, iz)) +
	     sqr(MRC_F3(fld, _RV1Y, ix, iy, iz)) +
	     sqr(MRC_F3(fld, _RV1Z, ix, iy, iz))) / MRC_F3(fld, _RR1, ix, iy, iz) +
      .5f * (sqr(MRC_F3(fld, _B1X, ix, iy, iz)) +
	     sqr(MRC_F3(fld, _B1Y, ix, iy, iz)) +
	     sqr(MRC_F3(fld, _B1Z, ix, iy, iz)));          
    

  } mrc_fld_foreach_end;
}

// ----------------------------------------------------------------------
// ggcm_mhd_ic_kh_descr

#define VAR(x) (void *)offsetof(struct ggcm_mhd_ic_kh, x)
static struct param ggcm_mhd_ic_kh_descr[] = {
  {"mpermi", VAR(mpermi), PARAM_FLOAT(1.)},
  {"pert", VAR(pert), PARAM_FLOAT(1e-7)},  
  {"r0", VAR(r0), PARAM_FLOAT(2.0)},
  {"r1", VAR(r1), PARAM_FLOAT(1.0)},
  {"v0", VAR(v0), PARAM_FLOAT(0.5)},
  {"v1", VAR(v1), PARAM_FLOAT(-0.5)},
  {"B0", VAR(B0), PARAM_FLOAT(0.0)}, 
  {"p0", VAR(p0), PARAM_FLOAT(2.5)},
  {"lambda", VAR(lambda), PARAM_FLOAT(4.0)},
  {},
};
#undef VAR

// ----------------------------------------------------------------------
// ggcm_mhd_ic_kh_ops

struct ggcm_mhd_ic_ops ggcm_mhd_ic_kh_ops = {
  .name        = "kh",
  .size        = sizeof(struct ggcm_mhd_ic_kh),
  .param_descr = ggcm_mhd_ic_kh_descr,
  .run         = ggcm_mhd_ic_kh_run,
};
