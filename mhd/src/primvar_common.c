
#include "ggcm_mhd_defs.h"
#include "ggcm_mhd_defs_extra.h"
#include "ggcm_mhd_private.h"

#include <mrc_profile.h>

#include <math.h>
#include <assert.h>

// ----------------------------------------------------------------------
// primvar_c

void
primvar_c(struct ggcm_mhd *mhd, int m_curr)
{
  static int PR;
  if (!PR) {
    PR = prof_register("primvar_c", 1., 0, 0);
  }
  prof_start(PR);

  int mhd_type;
  mrc_fld_get_param_int(mhd->fld, "mhd_type", &mhd_type);
  assert(mhd_type == MT_SEMI_CONSERVATIVE_GGCM ||
	 mhd_type == MT_SEMI_CONSERVATIVE);

  struct mrc_fld *f = mrc_fld_get_as(mhd->fld, FLD_TYPE);
  mrc_fld_data_t gamm = mhd->par.gamm;
  mrc_fld_data_t s = gamm - 1.f;

  mrc_fld_foreach(f, ix,iy,iz, 2, 2) {
    F3(f,_RR, ix,iy,iz) = F3(f, m_curr + _RR1, ix,iy,iz);
    mrc_fld_data_t rri = 1.f / F3(f, m_curr + _RR1, ix,iy,iz);
    F3(f,_VX, ix,iy,iz) = rri * F3(f, m_curr + _RV1X, ix,iy,iz);
    F3(f,_VY, ix,iy,iz) = rri * F3(f, m_curr + _RV1Y, ix,iy,iz);
    F3(f,_VZ, ix,iy,iz) = rri * F3(f, m_curr + _RV1Z, ix,iy,iz);
    mrc_fld_data_t rvv =
      F3(f,_VX, ix,iy,iz) * F3(f, m_curr + _RV1X, ix,iy,iz) +
      F3(f,_VY, ix,iy,iz) * F3(f, m_curr + _RV1Y, ix,iy,iz) +
      F3(f,_VZ, ix,iy,iz) * F3(f, m_curr + _RV1Z, ix,iy,iz);
    F3(f,_PP, ix,iy,iz) = s * (F3(f, m_curr + _UU1, ix,iy,iz) - .5f * rvv);
    mrc_fld_data_t cs2 = fmaxf(gamm * F3(f,_PP, ix,iy,iz) * rri, 0.f);
    F3(f,_CMSV, ix,iy,iz) = sqrtf(rvv * rri) + sqrtf(cs2);
  } mrc_fld_foreach_end;

  mrc_fld_put_as(f, mhd->fld);

  prof_stop(PR);
}

