#include "libggcm/ggcm_mhd_ic_private.h"
#include "ggcmtest.h"
#include <mrc_domain.h>
#include <mrc_crds.h>
#include <math.h>
#include <string.h>
#include <assert.h>

#define ggcm_mhd_cweno(obj) mrc_to_subobj(obj, struct mhd)

// ======================================================================
// ggcm_mhd_ic subclass "doubletearing"

struct ggcm_mhd_ic_doubletearing {
  float n0; // initial density 
};
// ----------------------------------------------------------------------
// ggcm_mhd_ic_doubletearing_run

static void
ggcm_mhd_ic_doubletearing_run(struct ggcm_mhd_ic *ic)
{
  struct ggcm_mhd_ic_doubletearing *sub = mrc_to_subobj(ic, struct ggcm_mhd_ic_doubletearing);
  struct ggcm_mhd *gmhd = ic->mhd;  
  struct mrc_fld *fld = gmhd->fld;
  struct mrc_crds *crds = mrc_domain_get_crds(gmhd->domain);  
  float xl[3], xh[3], L[3], r[3];
  mrc_crds_get_xl_xh(crds, xl, xh);
  for(int i=0; i<3; i++){
    L[i] = xh[i] - xl[i];
  }
  float gamma = gmhd->par.gamm;
  mrc_fld_foreach(fld, ix, iy, iz, 1, 1) {
    r[0] = MRC_CRD(crds, 0, ix);
    r[1] = MRC_CRD(crds, 1, iy);
    r[2] = MRC_CRD(crds, 2, iz);

    float kx = 2. * M_PI / L[0], ky = 2. * M_PI / L[1];
    MRC_F3(fld_psi, 0, ix,iy,iz) = ((mhd->pa) * L[0] / 
				    (4. * M_PI)) * (1. - cos(2*kx*r[0])) * sin(ky*r[1]);    
    
    MRC_F3(fld, _RR1, ix, iy, iz) = sub->n0+.5 * sqr(sin(kx * r[0]));
    MRC_F3(fld, _B1Y , ix, iy, iz) = cos(kx*r[0]);
    MRC_F3(fld, _B1X, ix,iy,iz) += 
      (MRC_F3(fld_psi,0, ix,iy+1,iz) - MRC_F3(fld_psi,0, ix,iy,iz)) /
      (MRC_CRD(crds,1,iy+1) - MRC_CRD(crds,1,iy) );    
    MRC_F3(fld, _B1Y, ix,iy,iz) += -
      (MRC_F3(fld_psi,0, ix+1,iy+1,iz) - MRC_F3(fld_psi,0, ix,iy+1,iz)) /
      (MRC_CRD(crds,0,ix+1) - MRC_CRD(crds,0, ix));    

    MRC_F3(fld, _UU1 , ix, iy, iz) = (sub->n0+.5 * sqr(sin(kx * r[0]))) /(gamma - 1.f) +
      .5f * (sqr(MRC_F3(fld, _RV1X, ix, iy, iz)) +
	     sqr(MRC_F3(fld, _RV1Y, ix, iy, iz)) +
	     sqr(MRC_F3(fld, _RV1Z, ix, iy, iz))) / MRC_F3(fld, _RR1, ix, iy, iz) +
      .5f * (sqr(MRC_F3(fld, _B1X, ix, iy, iz)) +
	     sqr(MRC_F3(fld, _B1Y, ix, iy, iz)) +
	     sqr(MRC_F3(fld, _B1Z, ix, iy, iz)));          
  } mrc_fld_foreach_end;
}

// ----------------------------------------------------------------------
// ggcm_mhd_ic_doubletearing_descr

#define VAR(x) (void *)offsetof(struct ggcm_mhd_ic_doubletearing, x)
static struct param ggcm_mhd_ic_doubletearing_descr[] = {
  {"n0", VAR(n0), PARAM_FLOAT(1.0)},
  {},
};
#undef VAR

// ----------------------------------------------------------------------
// ggcm_mhd_ic_doubletearing_ops

struct ggcm_mhd_ic_ops ggcm_mhd_ic_doubletearing_ops = {
  .name        = "doubletearing",
  .size        = sizeof(struct ggcm_mhd_ic_doubletearing),
  .param_descr = ggcm_mhd_ic_doubletearing_descr,
  .run         = ggcm_mhd_ic_doubletearing_run,
};
