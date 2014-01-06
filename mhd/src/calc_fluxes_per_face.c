#include "ggcm_mhd_step_cweno_private.h" 

#include "ggcm_mhd_private.h"
#include "ggcm_mhd_crds.h"
#include "ggcm_mhd_diag.h"
#include <mrc_domain.h>
#include <mrc_ddc.h>

#include <assert.h>
#include <stdio.h>
#include <math.h>

// ----------------------------------------------------------------------
// calc_fluxes_per_face
//
// this calculates fluxes on the face i using reconstructed variables (fld) that are
// given on the respective face
// (Ziegler 2004 section 3.1) 

void 
calc_fluxes_per_face(struct mrc_fld **_flux, struct ggcm_mhd *mhd, struct mrc_fld *_fld, int i)
{
  float gamma = mhd->par.gamm;

  struct mrc_fld *fld = mrc_fld_get_as(_fld, "float");
  struct mrc_fld *flux[3];
  for (int f = 0; f < 3; f++) {
    flux[f] = mrc_fld_get_as(_flux[f], "float");
  }

  mrc_fld_foreach(fld, ix, iy, iz, 1, 1) {

#if SEMICONSV
    float rhoi = 1.f / MRC_F3(fld, _RR1, ix,iy,iz);
      
    float pp = (gamma - 1.f) *
      (MRC_F3(fld, _UU1, ix,iy,iz) - .5f * rhoi * (sqr(MRC_F3(fld, _RV1X, ix,iy,iz)) +
						   sqr(MRC_F3(fld, _RV1Y, ix,iy,iz)) +
						   sqr(MRC_F3(fld, _RV1Z, ix,iy,iz))));

    // mass consv. 
    FLUX(flux, i, _RR1, ix,iy,iz) = MRC_F3(fld, _RV1X+i, ix,iy,iz);
    
    // momentum eq. 
    for (int j = 0; j < 3; j++) {
      FLUX(flux, j, _RV1X+i, ix,iy,iz) = 
	rhoi * MRC_F3(fld, _RV1X+j, ix,iy,iz) * MRC_F3(fld, _RV1X+i, ix,iy,iz) +
	((j == i) ? pp : 0.) ;
    }
    
    // energy eq. 
    FLUX(flux, i, _UU1, ix,iy,iz) =
      ((MRC_F3(fld, _UU1, ix,iy,iz) + pp) * MRC_F3(fld, _RV1X+i, ix,iy,iz)) * rhoi;

#else 
    float d_i = mhd->par.d_i;

    float rhoi = 1.f / MRC_F3(fld, _RR1, ix,iy,iz);
      
    float BB = (0.5f) * (sqr(MRC_F3(fld, _B1X, ix,iy,iz)) +
			 sqr(MRC_F3(fld, _B1Y, ix,iy,iz)) +
			 sqr(MRC_F3(fld, _B1Z, ix,iy,iz)));
    
    float mB = (MRC_F3(fld, _B1X, ix,iy,iz) * MRC_F3(fld, _RV1X, ix,iy,iz)) + 
      (MRC_F3(fld, _B1Y, ix,iy,iz) * MRC_F3(fld, _RV1Y, ix,iy,iz)) + 
      (MRC_F3(fld, _B1Z, ix,iy,iz) * MRC_F3(fld, _RV1Z, ix,iy,iz)) ; 

    
    float JB = -(MRC_F3(fld, _B1X, ix,iy,iz)  * MRC_F3(fld, _JX, ix,iy,iz))  
      -(MRC_F3(fld, _B1Y, ix,iy,iz) * MRC_F3(fld, _JY, ix,iy,iz))  
      -(MRC_F3(fld, _B1Z, ix,iy,iz) * MRC_F3(fld, _JZ, ix,iy,iz)) ; 

    float pp = (gamma - 1.f) *
      (MRC_F3(fld, _UU1, ix,iy,iz) - .5f * rhoi * (sqr(MRC_F3(fld, _RV1X, ix,iy,iz)) +
						 sqr(MRC_F3(fld, _RV1Y, ix,iy,iz)) +
						 sqr(MRC_F3(fld, _RV1Z, ix,iy,iz)))- 
       (.5f * (sqr(MRC_F3(fld, _B1X, ix,iy,iz)) +
	       sqr(MRC_F3(fld, _B1Y, ix,iy,iz)) +
	       sqr(MRC_F3(fld, _B1Z, ix,iy,iz)))));
    
    // mass consv. 
    FLUX(flux, i, _RR1, ix,iy,iz) = MRC_F3(fld, _RV1X+i, ix,iy,iz);
    
    // momentum eq. 
    for (int j = 0; j < 3; j++) {
      FLUX(flux, j, _RV1X+i, ix,iy,iz) = 
	rhoi * MRC_F3(fld, _RV1X+j, ix,iy,iz) * MRC_F3(fld, _RV1X+i, ix,iy,iz) +
	((j == i) ? pp : 0.) + 
	((j == i) ? BB : 0.) - (MRC_F3(fld, _B1X+i, ix,iy,iz) * MRC_F3(fld, _B1X+j, ix,iy,iz));
    }
    
    // energy eq. 
    FLUX(flux, i, _UU1, ix,iy,iz) =
      ( ((MRC_F3(fld, _UU1, ix,iy,iz) + pp + BB)*MRC_F3(fld, _RV1X+i, ix,iy,iz))-
	(mB * MRC_F3(fld, _B1X+i, ix,iy,iz)) + 
	(d_i * ( -0.5*MRC_F3(fld, _JX+i, ix,iy,iz)*BB - MRC_F3(fld, _B1X+i, ix,iy,iz)*JB)) ) * rhoi;
#endif 
  } mrc_fld_foreach_end;

  mrc_fld_put_as(fld, _fld);
  for (int f = 0; f < 3; f++) {
    mrc_fld_put_as(flux[f], _flux[f]);
  }
}
