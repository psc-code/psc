
#include "ggcm_mhd_ic_private.h"

#include "ggcm_mhd_private.h"
#include "ggcm_mhd_flds.h"
#include "ggcmtest.h"

#include <mrc_domain.h>
#include <mrc_crds.h>
#include <math.h>
#include <string.h>
#include <assert.h>

#define ggcm_mhd_cweno(obj) mrc_to_subobj(obj, struct mhd)

// ======================================================================
// ggcm_mhd_ic subclass "sod"

struct ggcm_mhd_ic_sod {
  float mpermi;  
  float initrad; // inital radius
  float pin; // initial inside  pressure
  float pout; // initial outside pressure
  float n0; // initial density 
  const char* pdim; 
};
// ----------------------------------------------------------------------
// ggcm_mhd_ic_sod_run

static void
ggcm_mhd_ic_sod_run(struct ggcm_mhd_ic *ic)
{
  struct ggcm_mhd_ic_sod *sub = mrc_to_subobj(ic, struct ggcm_mhd_ic_sod);
  struct ggcm_mhd *gmhd = ic->mhd;  
  struct mrc_f3 *f3 = ggcm_mhd_flds_get_mrc_f3(gmhd->flds_base);
  struct mrc_crds *crds = mrc_domain_get_crds(gmhd->domain);  
  float xl[3], xh[3], L[3], r[3];
  mrc_crds_get_xl_xh(crds, xl, xh);
  for(int i=0; i<3; i++){
    L[i] = xh[i] - xl[i];
  }
  float gamma = gmhd->par.gamm;
  mrc_f3_foreach(f3, ix, iy, iz, 1, 1) {
    r[0] = MRC_CRD(crds, 0, ix);
    r[1] = MRC_CRD(crds, 1, iy);
    r[2] = MRC_CRD(crds, 2, iz);
  
    if((strcmp(sub->pdim, "xy") || (strcmp(sub->pdim,"yx")))== 0){
	MRC_F3(f3, _RR1, ix, iy, iz) = sub->n0;
	if( sqrt((r[0]*r[0]) + (r[1]*r[1])) <= sub->initrad ){
	  MRC_F3(f3, _UU1 , ix, iy, iz) = sub->pin / (gamma - 1.f) +
	    .5f * (sqr(MRC_F3(f3, _RV1X, ix, iy, iz)) +
		   sqr(MRC_F3(f3, _RV1Y, ix, iy, iz)) +
		   sqr(MRC_F3(f3, _RV1Z, ix, iy, iz))) / MRC_F3(f3, _RR1, ix, iy, iz) +
	    .5f * (sqr(MRC_F3(f3, _B1X, ix, iy, iz)) +
		   sqr(MRC_F3(f3, _B1Y, ix, iy, iz)) +
		   sqr(MRC_F3(f3, _B1Z, ix, iy, iz)));      
	} else{	
	  MRC_F3(f3, _UU1 , ix, iy, iz) = sub->pout / (gamma - 1.f) +
	    .5f * (sqr(MRC_F3(f3, _RV1X, ix, iy, iz)) +
		   sqr(MRC_F3(f3, _RV1Y, ix, iy, iz)) +
		   sqr(MRC_F3(f3, _RV1Z, ix, iy, iz))) / MRC_F3(f3, _RR1, ix, iy, iz) +
	    .5f * (sqr(MRC_F3(f3, _B1X, ix, iy, iz)) +
		   sqr(MRC_F3(f3, _B1Y, ix, iy, iz)) +
		   sqr(MRC_F3(f3, _B1Z, ix, iy, iz)));      
	}
	MRC_F3(f3, _RV1X, ix, iy, iz) = 0.0;
	MRC_F3(f3, _RV1Y, ix, iy, iz) = 0.0;
	MRC_F3(f3, _RV1Z, ix, iy, iz) = 0.0;      
    } else if((strcmp(sub->pdim, "yz") || (strcmp(sub->pdim,"zy"))) == 0){
	  MRC_F3(f3, _RR1, ix, iy, iz) = sub->n0;
	  if( sqrt((r[1]*r[1]) + (r[2]*r[2])) <= sub->initrad ){	
	    MRC_F3(f3, _UU1 , ix, iy, iz) = sub->pin / (gamma - 1.f) +
	      .5f * (sqr(MRC_F3(f3, _RV1X, ix, iy, iz)) +
		     sqr(MRC_F3(f3, _RV1Y, ix, iy, iz)) +
		     sqr(MRC_F3(f3, _RV1Z, ix, iy, iz))) / MRC_F3(f3, _RR1, ix, iy, iz) +
	      .5f * (sqr(MRC_F3(f3, _B1X, ix, iy, iz)) +
		     sqr(MRC_F3(f3, _B1Y, ix, iy, iz)) +
		     sqr(MRC_F3(f3, _B1Z, ix, iy, iz)));      
	  } else{	
	    MRC_F3(f3, _UU1 , ix, iy, iz) = sub->pout / (gamma - 1.f) +
	      .5f * (sqr(MRC_F3(f3, _RV1X, ix, iy, iz)) +
		     sqr(MRC_F3(f3, _RV1Y, ix, iy, iz)) +
		     sqr(MRC_F3(f3, _RV1Z, ix, iy, iz))) / MRC_F3(f3, _RR1, ix, iy, iz) +
	      .5f * (sqr(MRC_F3(f3, _B1X, ix, iy, iz)) +
		     sqr(MRC_F3(f3, _B1Y, ix, iy, iz)) +
		     sqr(MRC_F3(f3, _B1Z, ix, iy, iz)));      
	  }
	  MRC_F3(f3, _RV1X, ix, iy, iz) = 0.0;
	  MRC_F3(f3, _RV1Y, ix, iy, iz) = 0.0;
	  MRC_F3(f3, _RV1Z, ix, iy, iz) = 0.0;      	  
      } else if((strcmp(sub->pdim, "xz") || (strcmp(sub->pdim,"zx"))) == 0){
	  MRC_F3(f3, _RR1, ix, iy, iz) = sub->n0;
	  if( sqrt((r[0]*r[0]) + (r[2]*r[2])) <= sub->initrad ){	
	    MRC_F3(f3, _UU1 , ix, iy, iz) = sub->pin / (gamma - 1.f) +
	      .5f * (sqr(MRC_F3(f3, _RV1X, ix, iy, iz)) +
		     sqr(MRC_F3(f3, _RV1Y, ix, iy, iz)) +
		     sqr(MRC_F3(f3, _RV1Z, ix, iy, iz))) / MRC_F3(f3, _RR1, ix, iy, iz) +
	      .5f * (sqr(MRC_F3(f3, _B1X, ix, iy, iz)) +
		     sqr(MRC_F3(f3, _B1Y, ix, iy, iz)) +
		     sqr(MRC_F3(f3, _B1Z, ix, iy, iz)));      
	  } else{	
	    MRC_F3(f3, _UU1 , ix, iy, iz) = sub->pout / (gamma - 1.f) +
	      .5f * (sqr(MRC_F3(f3, _RV1X, ix, iy, iz)) +
		     sqr(MRC_F3(f3, _RV1Y, ix, iy, iz)) +
		     sqr(MRC_F3(f3, _RV1Z, ix, iy, iz))) / MRC_F3(f3, _RR1, ix, iy, iz) +
	      .5f * (sqr(MRC_F3(f3, _B1X, ix, iy, iz)) +
		     sqr(MRC_F3(f3, _B1Y, ix, iy, iz)) +
		     sqr(MRC_F3(f3, _B1Z, ix, iy, iz)));      
	  }
	  MRC_F3(f3, _RV1X, ix, iy, iz) = 0.0;
	  MRC_F3(f3, _RV1Y, ix, iy, iz) = 0.0;
	  MRC_F3(f3, _RV1Z, ix, iy, iz) = 0.0;      	  
    } else {           
	  assert(0); /* unknown initial condition */
    }
  } mrc_f3_foreach_end;
}




// ----------------------------------------------------------------------
// ggcm_mhd_ic_sod_descr

#define VAR(x) (void *)offsetof(struct ggcm_mhd_ic_sod, x)
static struct param ggcm_mhd_ic_sod_descr[] = {
  {"mpermi", VAR(mpermi), PARAM_FLOAT(1.)},
  {"pin", VAR(pin), PARAM_FLOAT(10.0)},
  {"pout", VAR(pout), PARAM_FLOAT(0.1)},
  {"n0", VAR(n0), PARAM_FLOAT(1.0)},
  {"pdim", VAR(pdim), PARAM_STRING("xy")},  
  {},
};
#undef VAR

// ----------------------------------------------------------------------
// ggcm_mhd_ic_sod_ops

struct ggcm_mhd_ic_ops ggcm_mhd_ic_sod_ops = {
  .name        = "sod",
  .size        = sizeof(struct ggcm_mhd_ic_sod),
  .param_descr = ggcm_mhd_ic_sod_descr,
  .run         = ggcm_mhd_ic_sod_run,
};
