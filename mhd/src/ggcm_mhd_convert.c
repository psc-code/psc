
#include "ggcm_mhd_private.h"
#include "ggcm_mhd_defs.h"

#include <mrc_domain.h>

#include <string.h>

// ----------------------------------------------------------------------
// ggcm_mhd_convert_fc_from_primitive
//
// converts from primitive variables to fully-conservative in-place.
// No ghost points are set, the staggered B fields need to exist on all faces
// (that means one more than cell-centered dims)

void
ggcm_mhd_convert_fc_from_primitive(struct ggcm_mhd *mhd, struct mrc_fld *fld_base)
{
  float gamma_m1 = mhd->par.gamm - 1.;

  assert(strcmp(mrc_fld_type(fld_base), "mhd_fc_float") == 0);
  struct mrc_fld *fld = mrc_fld_get_as(fld_base, "mhd_fc_float");

  // don't go into ghost cells in invariant directions
  int gdims[3];
  mrc_domain_get_global_dims(mhd->domain, gdims);
  int dx = (gdims[0] == 1) ? 0 : 1;
  int dy = (gdims[1] == 1) ? 0 : 1;
  int dz = (gdims[2] == 1) ? 0 : 1;

  mrc_fld_foreach(fld, ix,iy,iz, 0, 0) {
    RV1X(fld, ix,iy,iz) = RR1(fld, ix,iy,iz) * V1X(fld, ix,iy,iz);
    RV1Y(fld, ix,iy,iz) = RR1(fld, ix,iy,iz) * V1Y(fld, ix,iy,iz);
    RV1Z(fld, ix,iy,iz) = RR1(fld, ix,iy,iz) * V1Z(fld, ix,iy,iz);
    UU1(fld, ix,iy,iz) = PP1(fld, ix,iy,iz) / gamma_m1
      + .5*(sqr(.5*(B1X(fld, ix,iy,iz) + B1X(fld, ix+dx,iy,iz))) + 
	    sqr(.5*(B1Y(fld, ix,iy,iz) + B1Y(fld, ix,iy+dy,iz))) +
	    sqr(.5*(B1Z(fld, ix,iy,iz) + B1Z(fld, ix,iy,iz+dz)))) 
      + .5*(sqr(RV1X(fld, ix,iy,iz)) +
	    sqr(RV1Y(fld, ix,iy,iz)) +
	    sqr(RV1Z(fld, ix,iy,iz))) / RR1(fld, ix,iy,iz);
  } mrc_fld_foreach_end;
  
  mrc_fld_put_as(fld, fld_base);
}

