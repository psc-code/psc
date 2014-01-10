
#include "ggcm_mhd_defs.h"
#include "ggcm_mhd_private.h"

#include <mrc_fld.h>
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

  struct mrc_fld *f = mrc_fld_get_as(mhd->fld, "float");

  mrc_fld_foreach(f, ix,iy,iz, 1, 2) {
    MRC_F3(f,_BX, ix,iy,iz) = .5f * (MRC_F3(f, m_curr + _B1X, ix,iy,iz) +
				     MRC_F3(f, m_curr + _B1X, ix-1,iy,iz));
    MRC_F3(f,_BY, ix,iy,iz) = .5f * (MRC_F3(f, m_curr + _B1Y, ix,iy,iz) +
				     MRC_F3(f, m_curr + _B1Y, ix,iy-1,iz));
    MRC_F3(f,_BZ, ix,iy,iz) = .5f * (MRC_F3(f, m_curr + _B1Z, ix,iy,iz) +
				     MRC_F3(f, m_curr + _B1Z, ix,iy,iz-1));
  } mrc_fld_foreach_end;

  mrc_fld_put_as(f, mhd->fld);

  prof_stop(PR);
}
