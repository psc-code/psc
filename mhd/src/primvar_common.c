
#include "ggcm_mhd_defs.h"
#include "ggcm_mhd_defs_extra.h"
#include "ggcm_mhd_private.h"

#include <mrc_profile.h>

#include <math.h>
#include <assert.h>

#include "pde/pde_defs.h"

// mhd options

#define OPT_EQN OPT_EQN_MHD_SCONS

#include "pde/pde_setup.c"
#include "pde/pde_mhd_setup.c"

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

  pde_setup(f);
  pde_mhd_setup(mhd);

  fld3d_t p_f;
  fld3d_setup(&p_f, f);

  for (int p = 0; p < mrc_fld_nr_patches(f); p++) {
    fld3d_get(&p_f, p);

    fld3d_foreach(ix,iy,iz, 2, 2) {
      F3S(p_f,_RR, ix,iy,iz) = F3S(p_f, m_curr + _RR1, ix,iy,iz);
      mrc_fld_data_t rri = 1.f / F3S(p_f, m_curr + _RR1, ix,iy,iz);
      F3S(p_f,_VX, ix,iy,iz) = rri * F3S(p_f, m_curr + _RV1X, ix,iy,iz);
      F3S(p_f,_VY, ix,iy,iz) = rri * F3S(p_f, m_curr + _RV1Y, ix,iy,iz);
      F3S(p_f,_VZ, ix,iy,iz) = rri * F3S(p_f, m_curr + _RV1Z, ix,iy,iz);
      mrc_fld_data_t rvv =
	F3S(p_f,_VX, ix,iy,iz) * F3S(p_f, m_curr + _RV1X, ix,iy,iz) +
	F3S(p_f,_VY, ix,iy,iz) * F3S(p_f, m_curr + _RV1Y, ix,iy,iz) +
	F3S(p_f,_VZ, ix,iy,iz) * F3S(p_f, m_curr + _RV1Z, ix,iy,iz);
      F3S(p_f,_PP, ix,iy,iz) = s * (F3S(p_f, m_curr + _UU1, ix,iy,iz) - .5f * rvv);
      mrc_fld_data_t cs2 = fmaxf(gamm * F3S(p_f,_PP, ix,iy,iz) * rri, 0.f);
      F3S(p_f,_CMSV, ix,iy,iz) = sqrtf(rvv * rri) + sqrtf(cs2);
    } fld3d_foreach_end;

    fld3d_put(&p_f, p);
  }   
  mrc_fld_put_as(f, mhd->fld);
    
  prof_stop(PR);
}
