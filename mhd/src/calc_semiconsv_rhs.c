#include "ggcm_mhd_step_cweno_private.h" 

#include "ggcm_mhd_private.h"
#include "ggcm_mhd_crds.h"
#include "ggcm_mhd_diag.h"
#include <mrc_domain.h>
#include <mrc_ddc.h>

// ----------------------------------------------------------------------
// calc_semiconsv_rhs
//
// calculates rhs for semi-conservative mhd

void 
calc_semiconsv_rhs(struct ggcm_mhd *mhd, struct mrc_fld *_rhs, struct mrc_fld *_flux[3])
{

  struct mrc_crds *crds = mrc_domain_get_crds(mhd->domain);
  struct mrc_fld *rhs = mrc_fld_get_as(_rhs, "mhd_fc_float");
  struct mrc_fld *fld = mrc_fld_get_as(mhd->fld, "float");
  
  //struct mrc_fld *_J_cc[3] = mrc_fld_get_as(mhd, _J_cc, "float");
  //struct mrc_fld *J_cc = mrc_fld_get_as(_J_cc, "float");
 
  struct mrc_fld *_E_cc = ggcm_mhd_get_fields(mhd, "E_cc", 3);
  struct mrc_fld *E_cc = mrc_fld_get_as(_E_cc, "float");
  struct mrc_fld *_J_cc = ggcm_mhd_get_fields(mhd, "J_cc", 3);
  //struct mrc_fld *J_cc = mrc_fld_get_as(_J_cc, "float");

  struct mrc_fld *flux[3];
  for (int f = 0; f < 3; f++) {
    flux[f] = mrc_fld_get_as(_flux[f], "float");
  }

 // initialize cell center Electric field structure      

  mrc_fld_foreach(E_cc, ix,iy,iz, 1, 1) {     
    MRC_F3(E_cc, 0, ix,iy,iz) = .25f*(- FLUX(flux, 1, _EZ,ix,iy,iz)    // G8 --> (Balsara & Spicer 1999 notation)
				      - FLUX(flux, 1, _EZ,ix,iy+1,iz)  // G8
				      + FLUX(flux, 2, _EY,ix,iy,iz)    // H7
				      + FLUX(flux, 2, _EY,ix,iy,iz+1));// H7    

    MRC_F3(E_cc, 1, ix,iy,iz) = .25f*(- FLUX(flux, 2, _EX,ix ,iy,iz)   // H6 
				      - FLUX(flux, 2, _EX,ix,iy,iz+1)  // H6
				      + FLUX(flux, 0, _EZ,ix ,iy,iz)   // F8
				      + FLUX(flux, 0, _EZ,ix+1,iy,iz));// F8

    MRC_F3(E_cc, 2, ix,iy,iz) = .25f*(- FLUX(flux, 0, _EY,ix,iy,iz)    // F7
				      - FLUX(flux, 0, _EY,ix+1,iy,iz)  // F7
				      + FLUX(flux, 1, _EX,ix,iy,iz)    // G6
				      + FLUX(flux, 1, _EX,ix,iy+1,iz));// G6  
  } mrc_fld_foreach_end;
 
  // calculate cell centered J
  ggcm_mhd_calc_currcc( mhd, mhd->fld, _B1X, _J_cc );
 

  // calculate neg divg 
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

  
  // add JdotE source term  
 struct mrc_fld *J_cc = mrc_fld_get_as(_J_cc, "float");

 mrc_fld_foreach(rhs, ix, iy, iz, 0, 0) {
    MRC_F3(rhs, _UU1, ix, iy, iz) += 
      MRC_F3(E_cc, 0, ix, iy, iz) * MRC_F3(J_cc, 0, ix, iy, iz) + 
      MRC_F3(E_cc, 1, ix, iy, iz) * MRC_F3(J_cc, 1, ix, iy, iz) + 
      MRC_F3(E_cc, 2, ix, iy, iz) * MRC_F3(J_cc, 2, ix, iy, iz) ;   
  } mrc_fld_foreach_end; 

  // add JxB source term
  mrc_fld_foreach(rhs, ix, iy, iz, 0, 0) {    
    MRC_F3(rhs, _RV1X, ix, iy, iz) +=  
      MRC_F3(J_cc, 1, ix, iy, iz) * 0.5*(MRC_F3(fld, _B1Z, ix,iy,iz)+MRC_F3(fld, _B1Z, ix,iy,iz+1)) -   
      MRC_F3(J_cc, 2, ix, iy, iz) * 0.5*(MRC_F3(fld, _B1Y, ix,iy,iz)+MRC_F3(fld, _B1Y, ix,iy+1,iz));
    MRC_F3(rhs, _RV1Y, ix, iy, iz) -= 
      MRC_F3(J_cc, 0, ix, iy, iz) * 0.5*(MRC_F3(fld, _B1Z, ix,iy,iz)+MRC_F3(fld, _B1Z, ix,iy,iz+1)) -   
      MRC_F3(J_cc, 2, ix, iy, iz) * 0.5*(MRC_F3(fld, _B1X, ix,iy,iz)+MRC_F3(fld, _B1X, ix+1,iy,iz))    ;
    MRC_F3(rhs, _RV1Z, ix, iy, iz) +=  
      MRC_F3(J_cc, 0, ix, iy, iz) * 0.5*(MRC_F3(fld, _B1Y, ix,iy,iz)+MRC_F3(fld, _B1Y, ix,iy+1,iz)) - 
      MRC_F3(J_cc, 1, ix, iy, iz) * 0.5*(MRC_F3(fld, _B1X, ix,iy,iz)+MRC_F3(fld, _B1X, ix+1,iy,iz));
  } mrc_fld_foreach_end; 

  for (int f = 0; f < 3; f++) {
    mrc_fld_put_as(flux[f], _flux[f]);
  }

  mrc_fld_destroy(fld);
  mrc_fld_destroy(J_cc);
  mrc_fld_destroy(E_cc);
  mrc_fld_put_as(rhs, _rhs);
}
