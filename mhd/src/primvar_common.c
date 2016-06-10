
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

  for (int p = 0; p < mrc_fld_nr_patches(f); p++) {
    mrc_fld_foreach(f, ix,iy,iz, 2, 2) {
      M3(f,_RR, ix,iy,iz, p) = M3(f, m_curr + _RR1, ix,iy,iz, p);
      mrc_fld_data_t rri = 1.f / M3(f, m_curr + _RR1, ix,iy,iz, p);
      M3(f,_VX, ix,iy,iz, p) = rri * M3(f, m_curr + _RV1X, ix,iy,iz, p);
      M3(f,_VY, ix,iy,iz, p) = rri * M3(f, m_curr + _RV1Y, ix,iy,iz, p);
      M3(f,_VZ, ix,iy,iz, p) = rri * M3(f, m_curr + _RV1Z, ix,iy,iz, p);
      mrc_fld_data_t rvv =
	M3(f,_VX, ix,iy,iz, p) * M3(f, m_curr + _RV1X, ix,iy,iz, p) +
	M3(f,_VY, ix,iy,iz, p) * M3(f, m_curr + _RV1Y, ix,iy,iz, p) +
	M3(f,_VZ, ix,iy,iz, p) * M3(f, m_curr + _RV1Z, ix,iy,iz, p);
      M3(f,_PP, ix,iy,iz, p) = s * (M3(f, m_curr + _UU1, ix,iy,iz, p) - .5f * rvv);
      mrc_fld_data_t cs2 = fmaxf(gamm * M3(f,_PP, ix,iy,iz, p) * rri, 0.f);
      M3(f,_CMSV, ix,iy,iz, p) = sqrtf(rvv * rri) + sqrtf(cs2);
    } mrc_fld_foreach_end;
  }   
  mrc_fld_put_as(f, mhd->fld);
    
  prof_stop(PR);
}
