#include "ggcm_mhd_step_cweno_private.h" 

#include "ggcm_mhd_private.h"
#include "ggcm_mhd_crds.h"
#include "ggcm_mhd_diag.h"
#include <mrc_domain.h>
#include <mrc_ddc.h>
#include <mrc_fld_as_float.h>

// ----------------------------------------------------------------------
// calc_ct_rhs 
// 
// calculate induction equation rhs with constrained transport method
// (Ziegler 2004 section 3.3) 

void
calc_ct_rhs(struct ggcm_mhd *mhd, struct mrc_fld *rhs, struct mrc_fld *flux[3])
{  
  // compute edge centered electric fields by interpolation (iv)  here tmp_fld are edge centered 
  // so that e.g.  MRC_F3(tmp_fld, 0, 0, 0, 0) is E_x 0,-1/2,-1/2  
  //         i.e.  MRC_F3(tmp_fld, 0, ix,iy,iz) is E_x ix,iy-1/2,iz-1/2   
  //           and MRC_F3(tmp_fld, 1, 0, 0, 0) is E_y -1/2,0,-1/2  etc etc 
                       

  struct mrc_ddc *ddc = mrc_domain_get_ddc(mhd->domain);
  mrc_ddc_set_param_int(ddc, "max_n_fields", 3);
  struct mrc_crds *crds = mrc_domain_get_crds(mhd->domain);
  struct mrc_fld *E_ec = ggcm_mhd_get_fields(mhd, "E_ec", 3);

  //initialize cell edge center Electric field structure      
  mrc_fld_foreach(E_ec, ix,iy,iz, 2, 2) { 
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
  
  assert(mrc_fld_nr_patches(rhs) == 1);
  int p = 0;

  mrc_fld_foreach(rhs, ix, iy,  iz, 2, 2) {
    BX_(rhs, ix, iy, iz, p) =
      (-((MRC_F3(E_ec, 2, ix, iy+1, iz) - MRC_F3(E_ec, 2, ix, iy, iz)) /
	 (.5f*( MRC_CRDY(crds, iy+1) - MRC_CRDY(crds, iy-1)))) 
       +((MRC_F3(E_ec, 1, ix, iy, iz+1) - MRC_F3(E_ec, 1, ix, iy, iz)) /
	 (.5f*( MRC_CRDZ(crds, iz+1) - MRC_CRDZ(crds, iz-1))))); 
    
    BY_(rhs, ix, iy, iz, p) =
      (-((MRC_F3(E_ec, 0, ix, iy, iz+1) - MRC_F3(E_ec, 0, ix, iy, iz)) /
	 (.5f*( MRC_CRD(crds, 2, iz+1) - MRC_CRD(crds, 2, iz-1))))
       +((MRC_F3(E_ec, 2, ix+1, iy, iz) - MRC_F3(E_ec, 2, ix, iy, iz)) /
	 (.5f*(MRC_CRD(crds, 0, ix+1) - MRC_CRD(crds, 0, ix-1))))); 
    
    BZ_(rhs, ix, iy, iz, p) =
      (-((MRC_F3( E_ec, 1, ix+1, iy, iz) - MRC_F3(E_ec, 1, ix, iy, iz)) /
	 (.5f*( MRC_CRD(crds, 0, ix+1) - MRC_CRD(crds, 0, ix-1))))  
       +((MRC_F3( E_ec, 0, ix, iy+1, iz) - MRC_F3(E_ec, 0, ix, iy, iz)) /
	 (.5f*( MRC_CRD(crds, 1, iy+1) - MRC_CRD(crds, 1, iy-1))))); 
  } mrc_fld_foreach_end;



#if SWBND
  int nnx = mhd->img[0], nny = mhd->img[1], nnz = mhd->img[2];
  if (mrc_domain_get_neighbor_rank(mhd->domain, (int [3]) { -1, 0, 0}) < 0) {
    for (int iz = -2; iz < nnz-2; iz++) {
      for (int iy = -2; iy < nny-2; iy++) {
	B1X(rhs,-2, iy, iz) = 0.0; 
	B1Y(rhs,-2, iy, iz) = 0.0;
	B1Z(rhs,-2, iy, iz) = 0.0;

	B1X(rhs,-1, iy, iz) = 0.0; 
	B1Y(rhs,-1, iy, iz) = 0.0;
	B1Z(rhs,-1, iy, iz) = 0.0;

	B1X(rhs, 0, iy, iz) = 0.0; 
	B1Y(rhs, 0, iy, iz) = 0.0;
	B1Z(rhs, 0, iy, iz) = 0.0;

	B1X(rhs, 1, iy, iz) = 0.0; 
	B1Y(rhs, 1, iy, iz) = 0.0;
	B1Z(rhs, 1, iy, iz) = 0.0;
      }
    }
  }
#endif
  mrc_fld_destroy(E_ec);
}
