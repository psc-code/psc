
#include "ggcm_mhd_defs.h"
#include "ggcm_mhd_private.h"
#include "ggcm_mhd_crds.h"

#include <mrc_domain.h>
#include <mrc_fld_as_double.h>

#include <math.h>
#include <string.h>

// ----------------------------------------------------------------------
// ggcm_mhd_calc_curcc
// FIXME: "clearncurr"which toggles dipole field subtraction is not implemented 
// default in input.defines is false, and that is what is done here. 

void
ggcm_mhd_calc_currcc(struct ggcm_mhd *mhd, struct mrc_fld *fld, int m,
		     struct mrc_fld *currcc)
{
  struct mrc_patch_info info;
  mrc_domain_get_local_patch_info(mhd->domain, 0, &info);
  //const char *dptype = ggcm_dipole_type(mhd->dipole);
  struct mrc_fld *tmp = mrc_fld_duplicate(currcc);
  float *bd4x = ggcm_mhd_crds_get_crd(mhd->crds, 0, BD4);
  float *bd4y = ggcm_mhd_crds_get_crd(mhd->crds, 1, BD4);
  float *bd4z = ggcm_mhd_crds_get_crd(mhd->crds, 2, BD4);

  struct mrc_fld *f = mrc_fld_get_as(fld, FLD_TYPE);
  struct mrc_fld *t = mrc_fld_get_as(tmp, FLD_TYPE);
  struct mrc_fld *c = mrc_fld_get_as(currcc, FLD_TYPE);
  
  mrc_fld_foreach(tmp,ix,iy,iz, 1, 1) {
    // compute current on edge first    
    F3(t,0,ix,iy,iz) = 
      (F3(f, m+2, ix,iy+1,iz) - F3(f, m+2, ix,iy,iz)) * bd4y[iy] - 
      (F3(f, m+1, ix,iy,iz+1) - F3(f, m+1, ix,iy,iz)) * bd4z[iz] ;     
    F3(t,1,ix,iy,iz) = 
      (F3(f, m  , ix,iy,iz+1) - F3(f, m  , ix,iy,iz)) * bd4z[iz] -
      (F3(f, m+2, ix+1,iy,iz) - F3(f, m+2, ix,iy,iz)) * bd4x[ix];         
    F3(t,2, ix,iy,iz) = 
      (F3(f,m+1 , ix+1,iy,iz) - F3(f, m+1, ix,iy,iz)) * bd4x[ix] - 
      (F3(f,m   , ix,iy+1,iz) - F3(f, m  , ix,iy,iz)) * bd4y[iy]; 
  } mrc_fld_foreach_end;
    
  mrc_fld_foreach(currcc, ix,iy,iz, 1, 1) {
    // average to the center 
    // FIXME: note, this originally used zmask, not ymask
    F3(c, 0, ix,iy,iz) = 0.25f *  
      (F3(t, 0, ix,iy,iz  ) + F3(t, 0, ix,iy-1,iz  ) + 
       F3(t, 0, ix,iy,iz-1) + F3(t, 0, ix,iy-1,iz-1));
    F3(c, 1, ix,iy,iz) = 0.25f *
      (F3(t, 1, ix,iy,iz  ) + F3(t, 1, ix-1,iy,iz  ) + 
       F3(t, 1, ix,iy,iz-1) + F3(t, 1, ix-1,iy,iz-1));
    F3(c, 2, ix,iy,iz) = 0.25f *
      (F3(t, 2, ix,iy  ,iz) + F3(t, 2, ix-1,iy,iz  ) + 
       F3(t, 2, ix,iy-1,iz) + F3(t, 2, ix-1,iy-1,iz));
  } mrc_fld_foreach_end;     

  mrc_fld_put_as(c, currcc);
  mrc_fld_put_as(t, tmp);
  mrc_fld_put_as(f, fld);
  mrc_fld_destroy(tmp); 
}

