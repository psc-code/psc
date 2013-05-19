
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
// ggcm_mhd_ic subclass "bw"

struct ggcm_mhd_ic_bw {
  float mpermi;
  const char* pdim; 
};
// ----------------------------------------------------------------------
// ggcm_mhd_ic_bw_run

static void
ggcm_mhd_ic_bw_run(struct ggcm_mhd_ic *ic)
{
  struct ggcm_mhd_ic_bw *sub = mrc_to_subobj(ic, struct ggcm_mhd_ic_bw);
  struct ggcm_mhd *gmhd = ic->mhd;  
  struct mrc_fld *fld = ggcm_mhd_flds_get_mrc_fld(gmhd->flds_base);
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
  
  if(strcmp(sub->pdim, "x") == 1){
    if(fabs(r[0]) < 0.5*L[0]){
      // Left                         
      MRC_F3(fld, _RR1, ix, iy, iz) = 1.0;
      MRC_F3(fld, _B1X , ix, iy, iz) = 0.75;
      MRC_F3(fld, _B1Y , ix, iy, iz) = 1.0;
      MRC_F3(fld, _B1Z , ix, iy, iz) = 0.0;
      MRC_F3(fld, _UU1 , ix, iy, iz) = 1.0 / (gamma - 1.f) +
	.5f * (sqr(MRC_F3(fld, _RV1X, ix, iy, iz)) +
	       sqr(MRC_F3(fld, _RV1Y, ix, iy, iz)) +
	       sqr(MRC_F3(fld, _RV1Z, ix, iy, iz))) / MRC_F3(fld, _RR1, ix, iy, iz) +
	.5f * (sqr(MRC_F3(fld, _B1X, ix, iy, iz)) +
	       sqr(MRC_F3(fld, _B1Y, ix, iy, iz)) +
	       sqr(MRC_F3(fld, _B1Z, ix, iy, iz)));      
    }else{
      // Right
      MRC_F3(fld, _RR1, ix, iy, iz) = 0.125;
      MRC_F3(fld, _B1X , ix, iy, iz) = 0.75;
      MRC_F3(fld, _B1Y , ix, iy, iz) = -1.0;
      MRC_F3(fld, _B1Z , ix, iy, iz) = 0.0;
      MRC_F3(fld, _UU1 , ix, iy, iz) = 0.1 / (gamma - 1.f) +
	.5f * (sqr(MRC_F3(fld, _RV1X, ix, iy, iz)) +
	       sqr(MRC_F3(fld, _RV1Y, ix, iy, iz)) +
	       sqr(MRC_F3(fld, _RV1Z, ix, iy, iz))) / MRC_F3(fld, _RR1, ix, iy, iz) +
	.5f * (sqr(MRC_F3(fld, _B1X, ix, iy, iz)) +
	       sqr(MRC_F3(fld, _B1Y, ix, iy, iz)) +
	       sqr(MRC_F3(fld, _B1Z, ix, iy, iz)));      
    }
    MRC_F3(fld, _RV1Z , ix, iy, iz) = 0.0;    
    MRC_F3(fld, _RV1X , ix, iy, iz) = 0.0;
    MRC_F3(fld, _RV1Y , ix, iy, iz) = 0.0;

  } else if(strcmp(sub->pdim, "y") == 1){
    if(fabs(r[1]) < 0.5*L[1]){
      // Left 
      MRC_F3(fld, _RR1, ix, iy, iz) = 1.0;
      MRC_F3(fld, _B1Y , ix, iy, iz) = 0.75;
      MRC_F3(fld, _B1X , ix, iy, iz) = 1.0;
      MRC_F3(fld, _B1Z , ix, iy, iz) = 0.0;
      MRC_F3(fld, _UU1 , ix, iy, iz) = 1.0 / (gamma - 1.f) +
	.5f * (sqr(MRC_F3(fld, _RV1X, ix, iy, iz)) +
	       sqr(MRC_F3(fld, _RV1Y, ix, iy, iz)) +
	       sqr(MRC_F3(fld, _RV1Z, ix, iy, iz))) / MRC_F3(fld, _RR1, ix, iy, iz) +
	.5f * (sqr(MRC_F3(fld, _B1X, ix, iy, iz)) +
	       sqr(MRC_F3(fld, _B1Y, ix, iy, iz)) +
	       sqr(MRC_F3(fld, _B1Z, ix, iy, iz)));      
    }else{
      // Right
      MRC_F3(fld, _RR1, ix, iy, iz) = 0.125;
      MRC_F3(fld, _B1Y , ix, iy, iz) = 0.75;
      MRC_F3(fld, _B1X , ix, iy, iz) = -1.0;
      MRC_F3(fld, _B1Z , ix, iy, iz) = 0.0;
      MRC_F3(fld, _UU1 , ix, iy, iz) = 0.1 / (gamma - 1.f) +
	.5f * (sqr(MRC_F3(fld, _RV1X, ix, iy, iz)) +
	       sqr(MRC_F3(fld, _RV1Y, ix, iy, iz)) +
	       sqr(MRC_F3(fld, _RV1Z, ix, iy, iz))) / MRC_F3(fld, _RR1, ix, iy, iz) +
	.5f * (sqr(MRC_F3(fld, _B1X, ix, iy, iz)) +
	       sqr(MRC_F3(fld, _B1Y, ix, iy, iz)) +
	       sqr(MRC_F3(fld, _B1Z, ix, iy, iz)));      
    }
    MRC_F3(fld, _RV1Z , ix, iy, iz) = 0.0;
    MRC_F3(fld, _RV1X , ix, iy, iz) = 0.0;
    MRC_F3(fld, _RV1Y , ix, iy, iz) = 0.0;
    
  } else if(strcmp(sub->pdim, "z") == 1){
    if(fabs(r[2]) < 0.5*L[2]){
      // Left 
      MRC_F3(fld, _RR1, ix, iy, iz) = 1.0;
      MRC_F3(fld, _RV1X , ix, iy, iz) = 0.0;
      MRC_F3(fld, _RV1Y , ix, iy, iz) = 0.0;	
      MRC_F3(fld, _B1Y , ix, iy, iz) = 0.0;
      MRC_F3(fld, _B1Z , ix, iy, iz) = 0.75;
      MRC_F3(fld, _B1X , ix, iy, iz) = 1.0;
      MRC_F3(fld, _UU1 , ix, iy, iz) = 1.0 / (gamma - 1.f) +
	.5f * (sqr(MRC_F3(fld, _RV1X, ix, iy, iz)) +
	       sqr(MRC_F3(fld, _RV1Y, ix, iy, iz)) +
	       sqr(MRC_F3(fld, _RV1Z, ix, iy, iz))) / MRC_F3(fld, _RR1, ix, iy, iz) +
	.5f * (sqr(MRC_F3(fld, _B1X, ix, iy, iz)) +
	       sqr(MRC_F3(fld, _B1Y, ix, iy, iz)) +
	       sqr(MRC_F3(fld, _B1Z, ix, iy, iz)));      
    }else{
      // Right
      MRC_F3(fld, _RR1, ix, iy, iz) = 0.125;
      MRC_F3(fld, _B1Y , ix, iy, iz) = 0.0;
      MRC_F3(fld, _B1X , ix, iy, iz) = -1.0;
      MRC_F3(fld, _B1Z , ix, iy, iz) = 0.75;
      MRC_F3(fld, _UU1 , ix, iy, iz) = 0.1 / (gamma - 1.f) +
	.5f * (sqr(MRC_F3(fld, _RV1X, ix, iy, iz)) +
	       sqr(MRC_F3(fld, _RV1Y, ix, iy, iz)) +
	       sqr(MRC_F3(fld, _RV1Z, ix, iy, iz))) / MRC_F3(fld, _RR1, ix, iy, iz) +
	.5f * (sqr(MRC_F3(fld, _B1X, ix, iy, iz)) +
	       sqr(MRC_F3(fld, _B1Y, ix, iy, iz)) +
	       sqr(MRC_F3(fld, _B1Z, ix, iy, iz)));      
    }
    MRC_F3(fld, _RV1X , ix, iy, iz) = 0.0;
    MRC_F3(fld, _RV1Y , ix, iy, iz) = 0.0;	
    MRC_F3(fld, _RV1Z , ix, iy, iz) = 0.0;       
  } else {           
    assert(0); /* unknown initial condition */
  }
  } mrc_fld_foreach_end;
}




// ----------------------------------------------------------------------
// ggcm_mhd_ic_bw_descr

#define VAR(x) (void *)offsetof(struct ggcm_mhd_ic_bw, x)
static struct param ggcm_mhd_ic_bw_descr[] = {
  {"mpermi", VAR(mpermi), PARAM_FLOAT(1.)},
  {"pdim", VAR(pdim), PARAM_STRING("x")},  
  {},
};
#undef VAR

// ----------------------------------------------------------------------
// ggcm_mhd_ic_bw_ops

struct ggcm_mhd_ic_ops ggcm_mhd_ic_bw_ops = {
  .name        = "bw",
  .size        = sizeof(struct ggcm_mhd_ic_bw),
  .param_descr = ggcm_mhd_ic_bw_descr,
  .run         = ggcm_mhd_ic_bw_run,
};
