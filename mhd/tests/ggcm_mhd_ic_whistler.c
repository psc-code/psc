
#include "ggcm_mhd_ic_private.h"

#include "ggcm_mhd_defs.h"
#include "ggcm_mhd_private.h"
#include "ggcm_mhd_flds.h"

#include <mrc_domain.h>
#include <mrc_crds.h>
#include <math.h>

#define ggcm_mhd_cweno(obj) mrc_to_subobj(obj, struct mhd)

// ======================================================================
// ggcm_mhd_ic subclass "whistler"

struct ggcm_mhd_ic_whistler {
  float mpermi;
  float Boz;
  float pert; 
  float eps; 
  float n0; 
  float lambda;   
};
// ----------------------------------------------------------------------
// ggcm_mhd_ic_whistler_run

static void
ggcm_mhd_ic_whistler_run(struct ggcm_mhd_ic *ic)
{
  struct ggcm_mhd_ic_whistler *sub = mrc_to_subobj(ic, struct ggcm_mhd_ic_whistler);
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

    // Xin etal (2005)      
    float kk= (sub->lambda * 2.*M_PI) / L[2] ;
    const int *dims = mrc_fld_dims(fld);
    int nz = dims[2];
    float vp= 2.*M_PI*sub->Boz*nz/((sub->mpermi)*(sub->n0)*L[2]); 
    MRC_F3(fld, _B1X, ix,iy,iz) = (sub->pert) * vp * sin( kk*r[2] ) ;       
    MRC_F3(fld, _B1Y, ix,iy,iz) = -(sub->pert) * vp * cos( kk*r[2] ) ;   
    MRC_F3(fld, _B1Z, ix,iy,iz) = sub->Boz ; 
    MRC_F3(fld, _RV1X, ix,iy,iz) = (sub->pert) * sub->n0 * sin( kk*r[2] ) ;      
    MRC_F3(fld, _RV1Y, ix,iy,iz) = -(sub->pert) * sub->n0 * cos( kk*r[2] ) ;
    MRC_F3(fld, _RR1, ix, iy, iz) = sub->n0 ;       
    MRC_F3(fld, _UU1 , ix, iy, iz) =  MRC_F3(fld, _RR1, ix, iy, iz) / (gamma -1.f) + 	
      .5f * (sqr(MRC_F3(fld, _RV1X, ix, iy, iz)) +
	     sqr(MRC_F3(fld, _RV1Y, ix, iy, iz)) +
	     sqr(MRC_F3(fld, _RV1Z, ix, iy, iz))) / MRC_F3(fld, _RR1, ix, iy, iz) +
      .5f * (sqr(MRC_F3(fld, _B1X, ix, iy, iz)) +
	     sqr(MRC_F3(fld, _B1Y, ix, iy, iz)) +
	     sqr(MRC_F3(fld, _B1Z, ix, iy, iz)));    
  } mrc_fld_foreach_end;
}

// ----------------------------------------------------------------------
// ggcm_mhd_ic_whistler_descr

#define VAR(x) (void *)offsetof(struct ggcm_mhd_ic_whistler, x)
static struct param ggcm_mhd_ic_whistler_descr[] = {
  {"mpermi", VAR(mpermi), PARAM_FLOAT(1.0)},
  {"pert", VAR(pert), PARAM_FLOAT(1e-5)},
  {"Boz", VAR(Boz), PARAM_FLOAT(1.0)},
  {"n0", VAR(n0), PARAM_FLOAT(25.)},
  {"lambda", VAR(lambda), PARAM_FLOAT(4)},  
  {},
};
#undef VAR

// ----------------------------------------------------------------------
// ggcm_mhd_ic_whistler_ops

struct ggcm_mhd_ic_ops ggcm_mhd_ic_whistler_ops = {
  .name        = "whistler",
  .size        = sizeof(struct ggcm_mhd_ic_whistler),
  .param_descr = ggcm_mhd_ic_whistler_descr,
  .run         = ggcm_mhd_ic_whistler_run,
};
