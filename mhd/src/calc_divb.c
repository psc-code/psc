
#include "ggcm_mhd_defs.h"
#include "ggcm_mhd_private.h"
#include "ggcm_mhd_crds.h"

#include <mrc_fld.h>
#include <mrc_domain.h>

#include <stdio.h>
#include <math.h>
#include <assert.h>

// ----------------------------------------------------------------------
// ggcm_mhd_calc_divb

void
ggcm_mhd_calc_divb(struct ggcm_mhd *mhd, struct mrc_fld *fld, struct mrc_fld *divb)
{
  struct mrc_patch_info info;
  mrc_domain_get_local_patch_info(fld->_domain, 0, &info);

  float *bd2x = ggcm_mhd_crds_get_crd(mhd->crds, 0, BD2);
  float *bd2y = ggcm_mhd_crds_get_crd(mhd->crds, 1, BD2);
  float *bd2z = ggcm_mhd_crds_get_crd(mhd->crds, 2, BD2);

  assert(fld->_data_type == MRC_NT_FLOAT);
  struct mrc_fld *f = mrc_fld_get_as(fld, mrc_fld_type(fld));
  struct mrc_fld *d = mrc_fld_get_as(divb, "float");

  float max = 0.;
  mrc_fld_foreach(divb, ix,iy,iz, 0, 0) {
    MRC_F3(divb,0, ix,iy,iz) = 
      (B1X(f, ix+1,iy,iz) - B1X(f, ix,iy,iz)) / bd2x[ix] +
      (B1Y(f, ix,iy+1,iz) - B1Y(f, ix,iy,iz)) / bd2y[iy] +
      (B1Z(f, ix,iy,iz+1) - B1Z(f, ix,iy,iz)) / bd2z[iz];
    MRC_F3(divb,0, ix,iy,iz) *= MRC_F3(f,_YMASK, ix,iy,iz);

    // the incoming solar wind won't match and hence divb != 0 here
    if (info.off[0] == 0 && ix <= 0)
      continue;

    max = fmaxf(max, fabsf(MRC_F3(divb,0, ix,iy,iz)));
  } mrc_fld_foreach_end;

  mrc_fld_put_as(f, fld);
  mrc_fld_put_as(d, divb);

  mprintf("max divb = %g\n", max);
}

