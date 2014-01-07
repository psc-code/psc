#include "ggcm_mhd_step_cweno_private.h" 

#include "ggcm_mhd_private.h"
#include "ggcm_mhd_crds.h"
#include "ggcm_mhd_diag.h"
#include <mrc_domain.h>
#include <mrc_ddc.h>

#define F3 MRC_F3 // FIXME

// ----------------------------------------------------------------------
// calc_ct_rhs 
// 
// calculate induction equation rhs with constrained transport method
// (Ziegler 2004 section 3.3) 

void
calc_ct_rhs(struct ggcm_mhd *mhd, struct mrc_fld *_rhs, struct mrc_fld *_flux[3])
{  
  // compute edge centered electric fields by interpolation (iv)  here tmp_fld are edge centered 
  // so that e.g.  MRC_F3(tmp_fld, 0, 0, 0, 0) is E_x 0,-1/2,-1/2  
  //         i.e.  MRC_F3(tmp_fld, 0, ix,iy,iz) is E_x ix,iy-1/2,iz-1/2   
  //           and MRC_F3(tmp_fld, 1, 0, 0, 0) is E_y -1/2,0,-1/2  etc etc 
                       

  struct mrc_crds *crds = mrc_domain_get_crds(mhd->domain);
  struct mrc_fld *_E_ec = ggcm_mhd_get_fields(mhd, "E_ec", 3);

  struct mrc_fld *E_ec = mrc_fld_get_as(_E_ec, "float");
  struct mrc_fld *flux[3];
  for (int f = 0; f < 3; f++) {
    flux[f] = mrc_fld_get_as(_flux[f], "float");
  }

  //initialize cell edge center Electric field structure      
  mrc_fld_foreach(E_ec, ix,iy,iz, 1, 1) { 
    MRC_F3(E_ec, 0, ix,iy,iz) = .25f*(- FLUX(flux, 1, _EZ, ix  ,iy  ,iz  )
				      - FLUX(flux, 1, _EZ, ix  ,iy  ,iz-1) 
				      + FLUX(flux, 2, _EY, ix  ,iy  ,iz  )
				      + FLUX(flux, 2, _EY, ix  ,iy-1,iz  ));    
    MRC_F3(E_ec, 1, ix,iy,iz) = .25f*(- FLUX(flux, 2, _EX, ix  ,iy  ,iz  )
				      - FLUX(flux, 2, _EX, ix-1,iy  ,iz  )
				      + FLUX(flux, 0, _EZ, ix  ,iy  ,iz  )
				      + FLUX(flux, 0, _EZ, ix  ,iy  ,iz-1));
    MRC_F3(E_ec, 2, ix,iy,iz) = .25f*(- FLUX(flux, 0, _EY, ix  ,iy  ,iz  )
				      - FLUX(flux, 0, _EY, ix  ,iy-1,iz  )
				      + FLUX(flux, 1, _EX, ix  ,iy  ,iz  )
				      + FLUX(flux, 1, _EX, ix-1,iy  ,iz  ));    
  } mrc_fld_foreach_end;
  
  for (int f = 0; f < 3; f++) {
    mrc_fld_put_as(flux[f], _flux[f]);
  }
  mrc_fld_put_as(E_ec, _E_ec);

  //  fill_ghost_fld(mhd, E_ec);

  struct mrc_fld *rhs = mrc_fld_get_as(_rhs, "mhd_fc_float");
  E_ec = mrc_fld_get_as(_E_ec, "float");

  mrc_fld_foreach(rhs, ix, iy,  iz, 1, 1) {
    B1X(rhs, ix, iy, iz) =  
      (-((MRC_F3(E_ec, 2, ix, iy+1, iz) - MRC_F3(E_ec, 2, ix, iy, iz)) /
	 (.5f*( MRC_CRDY(crds, iy+1) - MRC_CRDY(crds, iy-1)))) 
       +((MRC_F3(E_ec, 1, ix, iy, iz+1) - MRC_F3(E_ec, 1, ix, iy, iz)) /
	 (.5f*( MRC_CRDZ(crds, iz+1) - MRC_CRDZ(crds, iz-1))))); 
    
    B1Y(rhs, ix, iy, iz) = 
      (-((MRC_F3(E_ec, 0, ix, iy, iz+1) - MRC_F3(E_ec, 0, ix, iy, iz)) /
	 (.5f*( MRC_CRD(crds, 2, iz+1) - MRC_CRD(crds, 2, iz-1))))
       +((MRC_F3(E_ec, 2, ix+1, iy, iz) - MRC_F3(E_ec, 2, ix, iy, iz)) /
	 (.5f*(MRC_CRD(crds, 0, ix+1) - MRC_CRD(crds, 0, ix-1))))); 
    
    B1Z(rhs, ix, iy, iz) = 
      (-((MRC_F3( E_ec, 1, ix+1, iy, iz) - MRC_F3(E_ec, 1, ix, iy, iz)) /
	 (.5f*( MRC_CRD(crds, 0, ix+1) - MRC_CRD(crds, 0, ix-1))))  
       +((MRC_F3( E_ec, 0, ix, iy+1, iz) - MRC_F3(E_ec, 0, ix, iy, iz)) /
	 (.5f*( MRC_CRD(crds, 1, iy+1) - MRC_CRD(crds, 1, iy-1))))); 
  } mrc_fld_foreach_end;

  mrc_fld_put_as(E_ec, _E_ec);
  mrc_fld_put_as(rhs, _rhs);

  mrc_fld_destroy(_E_ec);
}
