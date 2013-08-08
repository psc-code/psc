
#include "ggcm_mhd_ic_private.h"

#include "ggcm_mhd_defs.h"
#include "ggcm_mhd_private.h"

#include <mrc_domain.h>
#include <mrc_fld.h>
#include <mrc_crds.h>
#include <math.h>
#include <string.h>
#include <assert.h>

// ======================================================================
// ggcm_mhd_ic subclass "wave_alfven"

struct ggcm_mhd_ic_wave_alfven {
  float pert; // pertubation amplitude
  float B0; // initial density 
  float n0; // n0 
  float p0; // initial pressure 
  const char* pdim; 
};
// ----------------------------------------------------------------------
// ggcm_mhd_ic_wave_alfven_run

static void
ggcm_mhd_ic_wave_alfven_run(struct ggcm_mhd_ic *ic)
{
  struct ggcm_mhd_ic_wave_alfven *sub = mrc_to_subobj(ic, struct ggcm_mhd_ic_wave_alfven);
  struct ggcm_mhd *gmhd = ic->mhd;  
  struct mrc_fld *fld = gmhd->fld;
  struct mrc_crds *crds = mrc_domain_get_crds(gmhd->domain);  
  float xl[3], xh[3], r[3];
  mrc_crds_get_param_float3(crds, "l", xl);
  mrc_crds_get_param_float3(crds, "h", xh);
  float gamma = gmhd->par.gamm;
  mrc_fld_foreach(fld, ix, iy, iz, 1, 1) {
    r[0] = MRC_CRD(crds, 0, ix);
    r[1] = MRC_CRD(crds, 1, iy);
    r[2] = MRC_CRD(crds, 2, iz);
  
    if(strcmp(sub->pdim, "x") == 0){
      float k = 2. * M_PI / (xh[0] - xl[0]);

      MRC_F3(fld, _RR1  , ix, iy, iz) = sub->n0;
      MRC_F3(fld, _RV1X , ix, iy, iz) = 0.0; 
      MRC_F3(fld, _RV1Y , ix, iy, iz) = sub->pert * sin(k * r[0]);
      MRC_F3(fld, _RV1Z , ix, iy, iz) = 0.0;
      MRC_F3(fld, _B1X  , ix, iy, iz) = sub->B0;
      MRC_F3(fld, _B1Y  , ix, iy, iz) = sub->pert * sin(k * r[0] );
      MRC_F3(fld, _B1Z  , ix, iy, iz) = 0.0;
      
      MRC_F3(fld, _UU1 , ix, iy, iz) = sub->p0 / (gamma - 1.f) +
	.5f * (sqr(MRC_F3(fld, _RV1X, ix, iy, iz)) +
	       sqr(MRC_F3(fld, _RV1Y, ix, iy, iz)) +
	       sqr(MRC_F3(fld, _RV1Z, ix, iy, iz))) / MRC_F3(fld, _RR1, ix, iy, iz) +
	.5f * (sqr(MRC_F3(fld, _B1X, ix, iy, iz)) +
	       sqr(MRC_F3(fld, _B1Y, ix, iy, iz)) +
	       sqr(MRC_F3(fld, _B1Z, ix, iy, iz))); 
      
    } else if(strcmp(sub->pdim, "y") == 0){
      MRC_F3(fld, _RR1, ix, iy, iz) = sub->n0;
      MRC_F3(fld, _RV1Y , ix, iy, iz) = 0.0; 
      MRC_F3(fld, _RV1X , ix, iy, iz) = sub->pert*sin( 4. * M_PI * r[1] );
      MRC_F3(fld, _RV1Z , ix, iy, iz) = 0.0;
      MRC_F3(fld, _B1Y , ix, iy, iz) = sub->B0;
      MRC_F3(fld, _B1X , ix, iy, iz) = sub->pert*sin( 4. * M_PI * r[1] );
      MRC_F3(fld, _B1Z , ix, iy, iz) = 0.0;
      
      MRC_F3(fld, _UU1 , ix, iy, iz) = sub->p0 / (gamma - 1.f) +
	.5f * (sqr(MRC_F3(fld, _RV1X, ix, iy, iz)) +
	       sqr(MRC_F3(fld, _RV1Y, ix, iy, iz)) +
	       sqr(MRC_F3(fld, _RV1Z, ix, iy, iz))) / MRC_F3(fld, _RR1, ix, iy, iz) +
	.5f * (sqr(MRC_F3(fld, _B1X, ix, iy, iz)) +
	       sqr(MRC_F3(fld, _B1Y, ix, iy, iz)) +
	       sqr(MRC_F3(fld, _B1Z, ix, iy, iz)));

    } else if(strcmp(sub->pdim, "z") == 0){
      MRC_F3(fld, _RR1, ix, iy, iz) = sub->n0;
      MRC_F3(fld, _RV1Y , ix, iy, iz) = 0.0; 
      MRC_F3(fld, _RV1X , ix, iy, iz) = sub->pert*sin( 4. * M_PI * r[2] );
      MRC_F3(fld, _RV1Z , ix, iy, iz) = 0.0;
      MRC_F3(fld, _B1Z , ix, iy, iz) = sub->B0;
      MRC_F3(fld, _B1X , ix, iy, iz) = sub->pert*sin( 4. * M_PI * r[2] );
      MRC_F3(fld, _B1Y , ix, iy, iz) = 0.0;
      
      MRC_F3(fld, _UU1 , ix, iy, iz) = sub->p0 / (gamma - 1.f) +
	.5f * (sqr(MRC_F3(fld, _RV1X, ix, iy, iz)) +
	       sqr(MRC_F3(fld, _RV1Y, ix, iy, iz)) +
	       sqr(MRC_F3(fld, _RV1Z, ix, iy, iz))) / MRC_F3(fld, _RR1, ix, iy, iz) +
	.5f * (sqr(MRC_F3(fld, _B1X, ix, iy, iz)) +
	       sqr(MRC_F3(fld, _B1Y, ix, iy, iz)) +
	       sqr(MRC_F3(fld, _B1Z, ix, iy, iz)));
      
    } else {           
	  assert(0); /* unknown initial condition */
    }
  } mrc_fld_foreach_end;
}

// ----------------------------------------------------------------------
// ggcm_mhd_ic_wave_alfven_descr

#define VAR(x) (void *)offsetof(struct ggcm_mhd_ic_wave_alfven, x)
static struct param ggcm_mhd_ic_wave_alfven_descr[] = {
  { "pert"           , VAR(pert)          , PARAM_FLOAT(1e-3)    },
  { "B0"             , VAR(B0)            , PARAM_FLOAT(1.0)     },
  { "n0"             , VAR(n0)            , PARAM_FLOAT(1.0)     },
  { "p0"             , VAR(p0)            , PARAM_FLOAT(1.0)     }, 
  { "pdim"           , VAR(pdim)          , PARAM_STRING("x")    },  
  {},
};
#undef VAR

// ----------------------------------------------------------------------
// ggcm_mhd_ic_wave_alfven_ops

struct ggcm_mhd_ic_ops ggcm_mhd_ic_wave_alfven_ops = {
  .name        = "wave_alfven",
  .size        = sizeof(struct ggcm_mhd_ic_wave_alfven),
  .param_descr = ggcm_mhd_ic_wave_alfven_descr,
  .run         = ggcm_mhd_ic_wave_alfven_run,
};
