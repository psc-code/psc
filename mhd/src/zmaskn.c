
#include "ggcm_mhd_defs.h"
#include "ggcm_mhd_private.h"

#include <mrc_fld_as_float.h>
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
  float va02i = 1.f / sqr(mhd->par.speedlimit / mhd->par.vvnorm);
  float eps   = 1e-15f;

  mrc_fld_foreach(f, ix,iy,iz, 2, 2) {
    float bb = 
      sqr(F3(f,_BX, ix,iy,iz)) + 
      sqr(F3(f,_BY, ix,iy,iz)) +
      sqr(F3(f,_BZ, ix,iy,iz));
    float rrm = fmaxf(eps, bb * va02i);
    F3(f, _ZMASK, ix,iy,iz) = F3(f, _YMASK, ix,iy,iz) * 
      fminf(1.f, F3(f, _RR, ix,iy,iz) / rrm);
  } mrc_fld_foreach_end;

  mrc_fld_put_as(f, mhd->fld);
  prof_stop(PR);
}

