#include "ggcm_mhd_step_cweno_private.h" 

#include "ggcm_mhd_private.h"
#include "ggcm_mhd_crds.h"
#include "ggcm_mhd_diag.h"
#include <mrc_domain.h>
#include <mrc_ddc.h>

// ----------------------------------------------------------------------
// calc_neg_divg
//
// calculates negative divergence 

void 
calc_neg_divg(struct ggcm_mhd *mhd, struct mrc_fld *_rhs, struct mrc_fld *_flux[3])
{
  struct mrc_crds *crds = mrc_domain_get_crds(mhd->domain);

#if SEMICONSV
  struct mrc_fld *rhs = mrc_fld_get_as(_rhs, "float");
#else 
  struct mrc_fld *rhs = mrc_fld_get_as(_rhs, "mhd_fc_float"); 
#endif 

  struct mrc_fld *flux[3];
  for (int f = 0; f < 3; f++) {
    flux[f] = mrc_fld_get_as(_flux[f], "float");
  }

  for (int m = 0; m <= _UU1; m++) {
    mrc_fld_foreach(rhs, ix, iy, iz, 0, 0) {
      int ind[3] = { ix, iy, iz };
      
      MRC_F3(rhs, m, ix, iy, iz) = 0.;
      for(int i=0; i<3; i++) {
	int dind[3] = {0, 0, 0};
	dind[i] = 1;      
	
	MRC_F3(rhs, m, ix, iy, iz) -=
	  (FLUX(flux, i, m, ix+dind[0],iy+dind[1],iz+dind[2]) - 
	   FLUX(flux, i, m, ix,iy,iz))
	  / (MRC_CRD(crds, i, ind[i]+1) - MRC_CRD(crds, i, ind[i]));
	assert(isfinite(MRC_F3(rhs, m, ix, iy, iz)));
      }
    } mrc_fld_foreach_end; 
  }

  for (int f = 0; f < 3; f++) {
    mrc_fld_put_as(flux[f], _flux[f]);
  }
  mrc_fld_put_as(rhs, _rhs);
}
