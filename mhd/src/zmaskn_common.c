
#include "ggcm_mhd_defs.h"
#include "ggcm_mhd_defs_extra.h"
#include "ggcm_mhd_private.h"

#include <mrc_profile.h>

#include <math.h>
#include <assert.h>

// ----------------------------------------------------------------------
// zmaskn_c

void
zmaskn_c(struct ggcm_mhd *mhd)
{
  static int PR;
  if (!PR) {
    PR = prof_register("zmaskn_c", 1., 0, 0);
  }
  prof_start(PR);

  struct mrc_fld *f = mrc_fld_get_as(mhd->fld, FLD_TYPE);
  float va02i = 1.f / sqr(mhd->par.speedlimit / mhd->vvnorm);
  float eps   = 1e-15f;

  for (int p = 0; p < mrc_fld_nr_patches(f); p++) {
    mrc_fld_foreach(f, ix,iy,iz, 2, 2) {
      float bb = 
	sqr(M3(f,_BX, ix,iy,iz, p)) + 
	sqr(M3(f,_BY, ix,iy,iz, p)) +
	sqr(M3(f,_BZ, ix,iy,iz, p));
      float rrm = fmaxf(eps, bb * va02i);
      M3(f, _ZMASK, ix,iy,iz, p) = M3(f, _YMASK, ix,iy,iz, p) * 
	fminf(1.f, M3(f, _RR, ix,iy,iz, p) / rrm);
    } mrc_fld_foreach_end;
  }
  mrc_fld_put_as(f, mhd->fld);
  prof_stop(PR);
}

