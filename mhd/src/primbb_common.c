
#include "ggcm_mhd_defs.h"
#include "ggcm_mhd_defs_extra.h"
#include "ggcm_mhd_private.h"

#include <mrc_profile.h>

#include <math.h>
#include <assert.h>

// ----------------------------------------------------------------------
// primbb_c

void
primbb_c(struct ggcm_mhd *mhd, int m_curr)
{
  static int PR;
  if (!PR) {
    PR = prof_register("primbb_c", 1., 0, 0);
  }
  prof_start(PR);

  int mhd_type;
  mrc_fld_get_param_int(mhd->fld, "mhd_type", &mhd_type);
  assert(mhd_type == MT_SEMI_CONSERVATIVE_GGCM);

  struct mrc_fld *f = mrc_fld_get_as(mhd->fld, FLD_TYPE);

  for (int p = 0; p < mrc_fld_nr_patches(f); p++) {
    mrc_fld_foreach(f, ix,iy,iz, 1, 2) {
      M3(f,_BX, ix,iy,iz, p) = .5f * (M3(f, m_curr + _B1X, ix,iy,iz, p) +
				      M3(f, m_curr + _B1X, ix-1,iy,iz, p));
      M3(f,_BY, ix,iy,iz, p) = .5f * (M3(f, m_curr + _B1Y, ix,iy,iz, p) +
				      M3(f, m_curr + _B1Y, ix,iy-1,iz, p));
      M3(f,_BZ, ix,iy,iz, p) = .5f * (M3(f, m_curr + _B1Z, ix,iy,iz, p) +
				      M3(f, m_curr + _B1Z, ix,iy,iz-1, p));
    } mrc_fld_foreach_end;
  }

  mrc_fld_put_as(f, mhd->fld);

  prof_stop(PR);
}

