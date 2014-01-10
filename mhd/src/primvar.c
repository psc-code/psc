
#include "ggcm_mhd_defs.h"
#include "ggcm_mhd_private.h"

#include <mrc_fld.h>
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

  struct mrc_fld *f = mrc_fld_get_as(mhd->fld, "float");
  float gamm = mhd->par.gamm;
  float s = gamm - 1.f;

  mrc_fld_foreach(f, ix,iy,iz, 2, 2) {
    MRC_F3(f,_RR, ix,iy,iz) = MRC_F3(f, m_curr + _RR1, ix,iy,iz);
    float rri = 1.f / MRC_F3(f, m_curr + _RR1, ix,iy,iz);
    MRC_F3(f,_VX, ix,iy,iz) = rri * MRC_F3(f, m_curr + _RV1X, ix,iy,iz);
    MRC_F3(f,_VY, ix,iy,iz) = rri * MRC_F3(f, m_curr + _RV1Y, ix,iy,iz);
    MRC_F3(f,_VZ, ix,iy,iz) = rri * MRC_F3(f, m_curr + _RV1Z, ix,iy,iz);
    float rvv =
      MRC_F3(f,_VX, ix,iy,iz) * MRC_F3(f, m_curr + _RV1X, ix,iy,iz) +
      MRC_F3(f,_VY, ix,iy,iz) * MRC_F3(f, m_curr + _RV1Y, ix,iy,iz) +
      MRC_F3(f,_VZ, ix,iy,iz) * MRC_F3(f, m_curr + _RV1Z, ix,iy,iz);
    MRC_F3(f,_PP, ix,iy,iz) = s * (MRC_F3(f, m_curr + _UU1, ix,iy,iz) - .5f * rvv);
    float cs2 = fmaxf(gamm * MRC_F3(f,_PP, ix,iy,iz) * rri, 0.f);
    MRC_F3(f,_CMSV, ix,iy,iz) = sqrtf(rvv * rri) + sqrtf(cs2);
  } mrc_fld_foreach_end;

  mrc_fld_put_as(f, mhd->fld);

  prof_stop(PR);
}

// ----------------------------------------------------------------------
// primvar1_c

void
primvar1_c(struct ggcm_mhd *mhd)
{
  return primvar_c(mhd, _RR1);
}

