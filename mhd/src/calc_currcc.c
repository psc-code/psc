
#include "ggcm_mhd_defs.h"
#include "ggcm_mhd_private.h"
#include "ggcm_mhd_crds.h"

#include <mrc_domain.h>
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

  struct mrc_fld *f = mrc_fld_get_as(fld, "float");
  struct mrc_fld *c = mrc_fld_get_as(currcc, "float");

  // average B first to cell centers 
  mrc_fld_foreach(tmp,ix,iy,iz, 1, 1) {
    MRC_F3(tmp,0,ix,iy,iz) = .5f*(MRC_F3(f, m, ix,iy,iz) + MRC_F3(f, m, ix+1,iy,iz));
    MRC_F3(tmp,1,ix,iy,iz) = .5f*(MRC_F3(f, m+1, ix,iy,iz) + MRC_F3(f, m+1, ix,iy+1,iz));
    MRC_F3(tmp,2,ix,iy,iz) = .5f*(MRC_F3(f, m+2, ix,iy,iz) + MRC_F3(f, m+2, ix,iy,iz+1));
  } mrc_fld_foreach_end;

  mrc_fld_foreach(tmp,ix,iy,iz, 1, 1) {
    // compute current      
    MRC_F3(c,0,ix,iy,iz) = 
      .5f*(MRC_F3(tmp, 2, ix,iy+1,iz) - MRC_F3(tmp, 2, ix,iy-1,iz)) * bd4y[iy] - 
      .5f*(MRC_F3(tmp, 1, ix,iy,iz+1) - MRC_F3(tmp, 1, ix,iy,iz-1)) * bd4z[iz] ;     
    MRC_F3(c,1,ix,iy,iz) = 
      .5f*(MRC_F3(tmp, 0, ix,iy,iz+1) - MRC_F3(tmp, 0, ix,iy,iz-1)) * bd4z[iz] -
      .5f*(MRC_F3(tmp, 2, ix+1,iy,iz) - MRC_F3(tmp, 2, ix-1,iy,iz)) * bd4x[ix];         
    MRC_F3(c,2, ix,iy,iz) = 
      .5f*(MRC_F3(tmp, 1, ix+1,iy,iz) - MRC_F3(tmp, 1, ix-1,iy,iz)) * bd4x[ix] - 
      .5f*(MRC_F3(tmp, 0, ix,iy+1,iz) - MRC_F3(tmp, 0, ix,iy-1,iz)) * bd4y[iy];     
  } mrc_fld_foreach_end;

  mrc_fld_put_as(c, currcc);
  mrc_fld_put_as(tmp, tmp);
  mrc_fld_put_as(f, fld);
  mrc_fld_destroy(tmp); 
}


