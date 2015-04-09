
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
  int gdims[3];
  mrc_domain_get_global_dims(fld->_domain, gdims);
  int dx = (gdims[0] > 1), dy = (gdims[1] > 1), dz = (gdims[2] > 1);

  struct mrc_patch_info info;
  mrc_domain_get_local_patch_info(mhd->domain, 0, &info);
  struct mrc_fld *tmp = mrc_fld_duplicate(currcc);

  struct mrc_fld *f = ggcm_mhd_fld_get_as(fld, FLD_TYPE, MT_SEMI_CONSERVATIVE_GGCM,
					  BX, BX + 3);
  struct mrc_fld *t = mrc_fld_get_as(tmp, FLD_TYPE);
  struct mrc_fld *c = mrc_fld_get_as(currcc, FLD_TYPE);
  
  for (int p = 0; p < mrc_fld_nr_patches(t); p++) {
    float *bd4x = ggcm_mhd_crds_get_crd_p(mhd->crds, 0, BD4, p);
    float *bd4y = ggcm_mhd_crds_get_crd_p(mhd->crds, 1, BD4, p);
    float *bd4z = ggcm_mhd_crds_get_crd_p(mhd->crds, 2, BD4, p);

    mrc_fld_foreach(t,ix,iy,iz, 1, 1) {
      // compute current on edge first
      M3(t,0,ix,iy,iz, p) =
	(M3(f, m+2, ix,iy+dy,iz, p) - M3(f, m+2, ix,iy,iz, p)) * bd4y[iy] -
	(M3(f, m+1, ix,iy,iz+dz, p) - M3(f, m+1, ix,iy,iz, p)) * bd4z[iz];
      M3(t,1,ix,iy,iz, p) =
	(M3(f, m  , ix,iy,iz+dz, p) - M3(f, m  , ix,iy,iz, p)) * bd4z[iz] -
	(M3(f, m+2, ix+dx,iy,iz, p) - M3(f, m+2, ix,iy,iz, p)) * bd4x[ix];
      M3(t,2, ix,iy,iz, p) =
	(M3(f,m+1 , ix+dx,iy,iz, p) - M3(f, m+1, ix,iy,iz, p)) * bd4x[ix] -
	(M3(f,m   , ix,iy+dy,iz, p) - M3(f, m  , ix,iy,iz, p)) * bd4y[iy];
    } mrc_fld_foreach_end;
  }

  for (int p = 0; p < mrc_fld_nr_patches(c); p++) {
    mrc_fld_foreach(c, ix,iy,iz, 1, 1) {
      // average to the center
      M3(c, 0, ix,iy,iz, p) = 0.25f *
	(M3(t, 0, ix,iy,iz   , p) + M3(t, 0, ix,iy-dy,iz   , p) +
	 M3(t, 0, ix,iy,iz-dz, p) + M3(t, 0, ix,iy-dy,iz-dz, p));
      M3(c, 1, ix,iy,iz, p) = 0.25f *
	(M3(t, 1, ix,iy,iz   , p) + M3(t, 1, ix-dx,iy,iz   , p) +
	 M3(t, 1, ix,iy,iz-dz, p) + M3(t, 1, ix-dx,iy,iz-dz, p));
      M3(c, 2, ix,iy,iz, p) = 0.25f *
	(M3(t, 2, ix,iy   ,iz, p) + M3(t, 2, ix-dx,iy   ,iz, p) +
	 M3(t, 2, ix,iy-dy,iz, p) + M3(t, 2, ix-dx,iy-dy,iz, p));
    } mrc_fld_foreach_end;
  }

  mrc_fld_put_as(c, currcc);
  mrc_fld_put_as(t, tmp);
  ggcm_mhd_fld_put_as(f, fld, 0, 0);
  mrc_fld_destroy(tmp); 
}

