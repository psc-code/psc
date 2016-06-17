
#include "ggcm_mhd_step_private.h"
#include "ggcm_mhd_private.h"
#include "ggcm_mhd_defs.h"
#include "ggcm_mhd_defs_extra.h"
#include "ggcm_mhd_crds.h"
#include "ggcm_mhd_diag_private.h"

#include <mrc_domain.h>
#include <mrc_profile.h>

#include <math.h>
#include <string.h>

#include "pde/pde_defs.h"

// mhd options

#define OPT_EQN OPT_EQN_MHD_SCONS

#include "pde/pde_setup.c"
#include "pde/pde_mhd_setup.c"
#include "pde/pde_mhd_primvar.c"
#include "pde/pde_mhd_primbb.c"
#include "pde/pde_mhd_zmaskn.c"
#include "pde/pde_mhd_badval_checks.c"

// TODO:
// - handle various resistivity models
// - handle limit2, limit3
// - handle lowmask

#define REPS (1.e-10f)

enum {
  LIMIT_NONE,
  LIMIT_1,
};

static mrc_fld_data_t s_mhd_time;

static float *s_fd1x, *s_fd1y, *s_fd1z;

#define FX1X(ix) PDE_CRDX_CC(ix)
#define FX1Y(iy) PDE_CRDY_CC(iy)
#define FX1Z(iz) PDE_CRDZ_CC(iz)

#define FX2X(ix) (sqr(FX1X(ix)))
#define FX2Y(iy) (sqr(FX1Y(iy)))
#define FX2Z(iz) (sqr(FX1Z(iz)))

// FD1 is really different if 'legacy_fd1' is used
#if 1
#define FD1X(ix) (s_fd1x[ix])
#define FD1Y(iy) (s_fd1y[iy])
#define FD1Z(iz) (s_fd1z[iz])
#else
#define FD1X(ix) PDE_INV_DX(ix)
#define FD1Y(iy) PDE_INV_DY(iy)
#define FD1Z(iz) PDE_INV_DZ(iz)
#endif

#define BD1X(ix) PDE_INV_DXF(ix+1)
#define BD1Y(iy) PDE_INV_DYF(iy+1)
#define BD1Z(iz) PDE_INV_DZF(iz+1)

#define BD2X(ix) PDE_DX(ix)
#define BD2Y(iy) PDE_DY(iy)
#define BD2Z(iz) PDE_DZ(iz)

#define BD3X(ix) PDE_INV_DX(ix)
#define BD3Y(iy) PDE_INV_DY(iy)
#define BD3Z(iz) PDE_INV_DZ(iz)

#define BD4X(ix) BD1X(ix)
#define BD4Y(iy) BD1Y(iy)
#define BD4Z(iz) BD1Z(iz)


// ======================================================================

static void
rmaskn_c(fld3d_t p_f)
{
  mrc_fld_data_t diffco = s_diffco;

  fld3d_foreach(ix,iy,iz, 2, 2) {
    F3S(p_f,_RMASK, ix,iy,iz) = 0.f;
    if (FX1X(ix) < s_diff_swbnd)
      continue;
    if (iy + s_patch.off[1] < s_diff_obnd)
      continue;
    if (iz + s_patch.off[2] < s_diff_obnd)
      continue;
    if (ix + s_patch.off[0] >= s_gdims[0] - s_diff_obnd)
      continue;
    if (iy + s_patch.off[1] >= s_gdims[1] - s_diff_obnd)
      continue;
    if (iz + s_patch.off[2] >= s_gdims[2] - s_diff_obnd)
      continue;
    F3S(p_f, _RMASK, ix,iy,iz) = diffco * F3S(p_f, _ZMASK, ix,iy,iz);
  } fld3d_foreach_end;
}

static void
vgflrr_c(fld3d_t p_f)
{
  fld3d_foreach(ix,iy,iz, 2, 2) {
    mrc_fld_data_t a = F3S(p_f,_RR, ix,iy,iz);
    F3S(p_f,_TMP1, ix,iy,iz) = a * F3S(p_f,_VX, ix,iy,iz);
    F3S(p_f,_TMP2, ix,iy,iz) = a * F3S(p_f,_VY, ix,iy,iz);
    F3S(p_f,_TMP3, ix,iy,iz) = a * F3S(p_f,_VZ, ix,iy,iz);
  } fld3d_foreach_end;
}

static void
vgflrvx_c(fld3d_t p_f)
{
  fld3d_foreach(ix,iy,iz, 2, 2) {
    mrc_fld_data_t a = F3S(p_f,_RR, ix,iy,iz) * F3S(p_f,_VX, ix,iy,iz);
    F3S(p_f,_TMP1, ix,iy,iz) = a * F3S(p_f,_VX, ix,iy,iz);
    F3S(p_f,_TMP2, ix,iy,iz) = a * F3S(p_f,_VY, ix,iy,iz);
    F3S(p_f,_TMP3, ix,iy,iz) = a * F3S(p_f,_VZ, ix,iy,iz);
  } fld3d_foreach_end;
}

static void
vgflrvy_c(fld3d_t p_f)
{
  fld3d_foreach(ix,iy,iz, 2, 2) {
    mrc_fld_data_t a = F3S(p_f,_RR, ix,iy,iz) * F3S(p_f,_VY, ix,iy,iz);
    F3S(p_f,_TMP1, ix,iy,iz) = a * F3S(p_f,_VX, ix,iy,iz);
    F3S(p_f,_TMP2, ix,iy,iz) = a * F3S(p_f,_VY, ix,iy,iz);
    F3S(p_f,_TMP3, ix,iy,iz) = a * F3S(p_f,_VZ, ix,iy,iz);
  } fld3d_foreach_end;
}

static void
vgflrvz_c(fld3d_t p_f)
{
  fld3d_foreach(ix,iy,iz, 2, 2) {
    mrc_fld_data_t a = F3S(p_f,_RR, ix,iy,iz) * F3S(p_f,_VZ, ix,iy,iz);
    F3S(p_f,_TMP1, ix,iy,iz) = a * F3S(p_f,_VX, ix,iy,iz);
    F3S(p_f,_TMP2, ix,iy,iz) = a * F3S(p_f,_VY, ix,iy,iz);
    F3S(p_f,_TMP3, ix,iy,iz) = a * F3S(p_f,_VZ, ix,iy,iz);
  } fld3d_foreach_end;
}

static void
vgfluu_c(fld3d_t p_f)
{
  mrc_fld_data_t s = s_gamma / (s_gamma - 1.f);
  fld3d_foreach(ix,iy,iz, 2, 2) {
    mrc_fld_data_t ep = s * F3S(p_f,_PP, ix,iy,iz) +
      .5f * F3S(p_f,_RR, ix,iy,iz) * (sqr(F3S(p_f,_VX, ix,iy,iz)) + 
				      sqr(F3S(p_f,_VY, ix,iy,iz)) + 
				      sqr(F3S(p_f,_VZ, ix,iy,iz)));
    F3S(p_f,_TMP1, ix,iy,iz) = ep * F3S(p_f,_VX, ix,iy,iz);
    F3S(p_f,_TMP2, ix,iy,iz) = ep * F3S(p_f,_VY, ix,iy,iz);
    F3S(p_f,_TMP3, ix,iy,iz) = ep * F3S(p_f,_VZ, ix,iy,iz);
  } fld3d_foreach_end;
}

static void
fluxl_c(fld3d_t p_f, int m)
{
  fld3d_foreach(ix,iy,iz, 1, 0) {
    mrc_fld_data_t aa = F3S(p_f,m, ix,iy,iz);
    mrc_fld_data_t cmsv = F3S(p_f,_CMSV, ix,iy,iz);
    F3S(p_f,_FLX, ix,iy,iz) =
      .5f * ((F3S(p_f,_TMP1, ix  ,iy,iz) + F3S(p_f,_TMP1, ix+1,iy,iz)) -
	     .5f * (F3S(p_f,_CMSV, ix+1,iy,iz) + cmsv) * (F3S(p_f,m, ix+1,iy,iz) - aa));
    F3S(p_f,_FLY, ix,iy,iz) =
      .5f * ((F3S(p_f,_TMP2, ix,iy  ,iz) + F3S(p_f,_TMP2, ix,iy+1,iz)) -
	     .5f * (F3S(p_f,_CMSV, ix,iy+1,iz) + cmsv) * (F3S(p_f,m, ix,iy+1,iz) - aa));
    F3S(p_f,_FLZ, ix,iy,iz) =
      .5f * ((F3S(p_f,_TMP3, ix,iy,iz  ) + F3S(p_f,_TMP3, ix,iy,iz+1)) -
	     .5f * (F3S(p_f,_CMSV, ix,iy,iz+1) + cmsv) * (F3S(p_f,m, ix,iy,iz+1) - aa));
  } fld3d_foreach_end;
}

static void
fluxb_c(fld3d_t p_f, int m)
{
  mrc_fld_data_t s1 = 1.f/12.f;
  mrc_fld_data_t s7 = 7.f * s1;

  fld3d_foreach(ix,iy,iz, 1, 0) {
    mrc_fld_data_t fhx = (s7 * (F3S(p_f, _TMP1, ix  ,iy,iz) + F3S(p_f, _TMP1, ix+1,iy,iz)) -
		 s1 * (F3S(p_f, _TMP1, ix-1,iy,iz) + F3S(p_f, _TMP1, ix+2,iy,iz)));
    mrc_fld_data_t fhy = (s7 * (F3S(p_f, _TMP2, ix,iy  ,iz) + F3S(p_f, _TMP2, ix,iy+1,iz)) -
		 s1 * (F3S(p_f, _TMP2, ix,iy-1,iz) + F3S(p_f, _TMP2, ix,iy+2,iz)));
    mrc_fld_data_t fhz = (s7 * (F3S(p_f, _TMP3, ix,iy,iz  ) + F3S(p_f, _TMP3, ix,iy,iz+1)) -
		 s1 * (F3S(p_f, _TMP3, ix,iy,iz-1) + F3S(p_f, _TMP3, ix,iy,iz+2)));

    mrc_fld_data_t aa = F3S(p_f,m, ix,iy,iz);
    mrc_fld_data_t cmsv = F3S(p_f,_CMSV, ix,iy,iz);
    mrc_fld_data_t flx =
      .5f * ((F3S(p_f,_TMP1, ix  ,iy,iz) + F3S(p_f,_TMP1, ix+1,iy,iz)) -
	     .5f * (F3S(p_f,_CMSV, ix+1,iy,iz) + cmsv) * (F3S(p_f,m, ix+1,iy,iz) - aa));
    mrc_fld_data_t fly =
      .5f * ((F3S(p_f,_TMP2, ix,iy  ,iz) + F3S(p_f,_TMP2, ix,iy+1,iz)) -
	     .5f * (F3S(p_f,_CMSV, ix,iy+1,iz) + cmsv) * (F3S(p_f,m, ix,iy+1,iz) - aa));
    mrc_fld_data_t flz = 
      .5f * ((F3S(p_f,_TMP3, ix,iy,iz  ) + F3S(p_f,_TMP3, ix,iy,iz+1)) -
	     .5f * (F3S(p_f,_CMSV, ix,iy,iz+1) + cmsv) * (F3S(p_f,m, ix,iy,iz+1) - aa));

    mrc_fld_data_t cx = F3S(p_f, _CX, ix,iy,iz);
    F3S(p_f, _FLX, ix,iy,iz) = cx * flx + (1.f - cx) * fhx;
    mrc_fld_data_t cy = F3S(p_f, _CY, ix,iy,iz);
    F3S(p_f, _FLY, ix,iy,iz) = cy * fly + (1.f - cy) * fhy;
    mrc_fld_data_t cz = F3S(p_f, _CZ, ix,iy,iz);
    F3S(p_f, _FLZ, ix,iy,iz) = cz * flz + (1.f - cz) * fhz;
  } fld3d_foreach_end;
}

static void
pushn_c(fld3d_t p_f, int ma, int mc, mrc_fld_data_t dt)
{
  fld3d_foreach(ix,iy,iz, 0, 0) {
    mrc_fld_data_t s = dt * F3S(p_f,_YMASK, ix,iy,iz);
    F3S(p_f,mc, ix,iy,iz) = F3S(p_f,ma, ix,iy,iz)
      - s * (FD1X(ix) * (F3S(p_f,_FLX, ix,iy,iz) - F3S(p_f,_FLX, ix-1,iy,iz)) +
	     FD1Y(iy) * (F3S(p_f,_FLY, ix,iy,iz) - F3S(p_f,_FLY, ix,iy-1,iz)) +
	     FD1Z(iz) * (F3S(p_f,_FLZ, ix,iy,iz) - F3S(p_f,_FLZ, ix,iy,iz-1)));
  } fld3d_foreach_end;
}

static void
pushpp_c(fld3d_t p_f, mrc_fld_data_t dt, int m)
{
  mrc_fld_data_t dth = -.5f * dt;
  fld3d_foreach(ix,iy,iz, 0, 0) {
    mrc_fld_data_t fpx = FD1X(ix) * (F3S(p_f, _PP, ix+1,iy,iz) - F3S(p_f, _PP, ix-1,iy,iz));
    mrc_fld_data_t fpy = FD1Y(iy) * (F3S(p_f, _PP, ix,iy+1,iz) - F3S(p_f, _PP, ix,iy-1,iz));
    mrc_fld_data_t fpz = FD1Z(iz) * (F3S(p_f, _PP, ix,iy,iz+1) - F3S(p_f, _PP, ix,iy,iz-1));
    mrc_fld_data_t z = dth * F3S(p_f,_ZMASK, ix,iy,iz);
    F3S(p_f, m + _RV1X, ix,iy,iz) += z * fpx;
    F3S(p_f, m + _RV1Y, ix,iy,iz) += z * fpy;
    F3S(p_f, m + _RV1Z, ix,iy,iz) += z * fpz;
  } fld3d_foreach_end;
}

static void
vgrs(fld3d_t p_f, int m, mrc_fld_data_t s)
{
  fld3d_foreach(ix,iy,iz, 2, 2) {
    F3S(p_f, m, ix,iy,iz) = s;
  } fld3d_foreach_end;
}

static void
vgrv(fld3d_t p_f, int m_to, int m_from)
{
  fld3d_foreach(ix,iy,iz, 2, 2) {
    F3S(p_f, m_to, ix,iy,iz) = F3S(p_f, m_from, ix,iy,iz);
  } fld3d_foreach_end;
}

static inline void
limit1a(fld3d_t p_f, int m, int ix, int iy, int iz, int IX, int IY, int IZ, int C)
{
  const mrc_fld_data_t reps = 0.003;
  const mrc_fld_data_t seps = -0.001;
  const mrc_fld_data_t teps = 1.e-25;

  // Harten/Zwas type switch
  mrc_fld_data_t aa = F3S(p_f, m, ix,iy,iz);
  mrc_fld_data_t a1 = F3S(p_f, m, ix+IX,iy+IY,iz+IZ);
  mrc_fld_data_t a2 = F3S(p_f, m, ix-IX,iy-IY,iz-IZ);
  mrc_fld_data_t d1 = aa - a2;
  mrc_fld_data_t d2 = a1 - aa;
  mrc_fld_data_t s1 = mrc_fld_abs(d1);
  mrc_fld_data_t s2 = mrc_fld_abs(d2);
  mrc_fld_data_t f1 = mrc_fld_abs(a1) + mrc_fld_abs(a2) + mrc_fld_abs(aa);
  mrc_fld_data_t s5 = s1 + s2 + reps*f1 + teps;
  mrc_fld_data_t r3 = mrc_fld_abs(s1 - s2) / s5; // edge condition
  mrc_fld_data_t f2 = seps * f1 * f1;
  if (d1 * d2 < f2) {
    r3 = 1.f;
  }
  r3 = r3 * r3;
  r3 = r3 * r3;
  r3 = mrc_fld_min(2.f * r3, 1.);
  F3S(p_f, C, ix   ,iy   ,iz   ) = mrc_fld_max(F3S(p_f, C, ix   ,iy   ,iz   ), r3);
  F3S(p_f, C, ix-IX,iy-IY,iz-IZ) = mrc_fld_max(F3S(p_f, C, ix-IX,iy-IY,iz-IZ), r3);
}

static void
limit1_c(fld3d_t p_f, int m, int C)
{
  if (s_mhd_time < s_timelo) {
    vgrs(p_f, C + 0, 1.f);
    vgrs(p_f, C + 1, 1.f);
    vgrs(p_f, C + 2, 1.f);
    return;
  }

  fld3d_foreach(ix,iy,iz, 1, 1) {
/* .if (limit_aspect_low) then */
/* .call lowmask(0,0,0,tl1) */
    limit1a(p_f, m, ix,iy,iz, 1,0,0, C + 0);
    limit1a(p_f, m, ix,iy,iz, 0,1,0, C + 1);
    limit1a(p_f, m, ix,iy,iz, 0,0,1, C + 2);
  } fld3d_foreach_end;
}

static void
vgfl_c(fld3d_t p_f, int m)
{
  switch (m) {
  case _RR1:  return vgflrr_c(p_f);
  case _RV1X: return vgflrvx_c(p_f);
  case _RV1Y: return vgflrvy_c(p_f);
  case _RV1Z: return vgflrvz_c(p_f);
  case _UU1:  return vgfluu_c(p_f);
  default: assert(0);
  }
}

static void
pushfv_c(fld3d_t p_f, int m, mrc_fld_data_t dt,
	 int m_prev, int m_curr, int m_next, int limit)
{
  vgfl_c(p_f, m);
  if (limit == LIMIT_NONE) {
    fluxl_c(p_f, m_curr + m);
  } else {
    vgrv(p_f, _CX, _BX); vgrv(p_f, _CY, _BY); vgrv(p_f, _CZ, _BZ);
    limit1_c(p_f, m_curr + m, _CX);
    fluxb_c(p_f, m_curr + m);
  }

  pushn_c(p_f, m_prev + m, m_next + m, dt);
}

// ----------------------------------------------------------------------
// curr_c
//
// edge centered current density

static void
curr_c(fld3d_t p_f, int m_j, int m_curr)
{
  fld3d_foreach(ix,iy,iz, 2, 1) {
    F3S(p_f, m_j + 0, ix,iy,iz) =
      (F3S(p_f, m_curr + _B1Z, ix,iy+1,iz) - F3S(p_f, m_curr + _B1Z, ix,iy,iz)) * BD4Y(iy) -
      (F3S(p_f, m_curr + _B1Y, ix,iy,iz+1) - F3S(p_f, m_curr + _B1Y, ix,iy,iz)) * BD4Z(iz);
    F3S(p_f, m_j + 1, ix,iy,iz) =
      (F3S(p_f, m_curr + _B1X, ix,iy,iz+1) - F3S(p_f, m_curr + _B1X, ix,iy,iz)) * BD4Z(iz) -
      (F3S(p_f, m_curr + _B1Z, ix+1,iy,iz) - F3S(p_f, m_curr + _B1Z, ix,iy,iz)) * BD4X(ix);
    F3S(p_f, m_j + 2, ix,iy,iz) =
      (F3S(p_f, m_curr + _B1Y, ix+1,iy,iz) - F3S(p_f, m_curr + _B1Y, ix,iy,iz)) * BD4X(ix) -
      (F3S(p_f, m_curr + _B1X, ix,iy+1,iz) - F3S(p_f, m_curr + _B1X, ix,iy,iz)) * BD4Y(iy);
  } fld3d_foreach_end;
}

// ----------------------------------------------------------------------
// currbb_c
//
// cell-averaged B

static void
currbb_c(fld3d_t p_f, int m, int m_curr)
{
  fld3d_foreach(ix,iy,iz, 1, 1) {
    F3S(p_f, m+0, ix,iy,iz) = .5f * (F3S(p_f, m_curr + _B1X, ix  ,iy,iz) +
				     F3S(p_f, m_curr + _B1X, ix-1,iy,iz));
    F3S(p_f, m+1, ix,iy,iz) = .5f * (F3S(p_f, m_curr + _B1Y, ix,iy  ,iz) +
				     F3S(p_f, m_curr + _B1Y, ix,iy-1,iz));
    F3S(p_f, m+2, ix,iy,iz) = .5f * (F3S(p_f, m_curr + _B1Z, ix,iy,iz  ) +
				     F3S(p_f, m_curr + _B1Z, ix,iy,iz-1));
  } fld3d_foreach_end;
}

// ----------------------------------------------------------------------
// curbc_c
//
// cell-centered j

static void
curbc_c(fld3d_t p_f, int m_curr)
{ 
  enum { _TX = _TMP1, _TY = _TMP2, _TZ = _TMP3 };

  curr_c(p_f, _TX, m_curr);

  // j averaged to cell-centered
  fld3d_foreach(ix,iy,iz, 1, 1) {
    mrc_fld_data_t s = .25f * F3S(p_f, _ZMASK, ix, iy, iz);
    F3S(p_f, _CURRX, ix,iy,iz) = s * (F3S(p_f, _TX, ix,iy  ,iz  ) + F3S(p_f, _TX, ix,iy-1,iz  ) +
				      F3S(p_f, _TX, ix,iy  ,iz-1) + F3S(p_f, _TX, ix,iy-1,iz-1));
    F3S(p_f, _CURRY, ix,iy,iz) = s * (F3S(p_f, _TY, ix  ,iy,iz  ) + F3S(p_f, _TY, ix-1,iy,iz  ) +
				      F3S(p_f, _TY, ix  ,iy,iz-1) + F3S(p_f, _TY, ix-1,iy,iz-1));
    F3S(p_f, _CURRZ, ix,iy,iz) = s * (F3S(p_f, _TZ, ix  ,iy  ,iz) + F3S(p_f, _TZ, ix-1,iy  ,iz) +
				      F3S(p_f, _TZ, ix  ,iy-1,iz) + F3S(p_f, _TZ, ix-1,iy-1,iz));
  } fld3d_foreach_end;
}

static void
push_ej_c(fld3d_t p_f, mrc_fld_data_t dt, int m_curr, int m_next)
{
  enum { XJX = _BX, XJY = _BY, XJZ = _BZ };
  enum { BX = _TMP1, BY = _TMP2, BZ = _TMP3 };

  curr_c(p_f, XJX, m_curr);
  currbb_c(p_f, BX, m_curr);
	
  mrc_fld_data_t s1 = .25f * dt;
  fld3d_foreach(ix,iy,iz, 0, 0) {
    mrc_fld_data_t z = F3S(p_f,_ZMASK, ix,iy,iz);
    mrc_fld_data_t s2 = s1 * z;
    mrc_fld_data_t cx = (F3S(p_f, XJX, ix  ,iy  ,iz  ) +
		F3S(p_f, XJX, ix  ,iy-1,iz  ) +
		F3S(p_f, XJX, ix  ,iy  ,iz-1) +
		F3S(p_f, XJX, ix  ,iy-1,iz-1));
    mrc_fld_data_t cy = (F3S(p_f, XJY, ix  ,iy  ,iz  ) +
		F3S(p_f, XJY, ix-1,iy  ,iz  ) +
		F3S(p_f, XJY, ix  ,iy  ,iz-1) +
		F3S(p_f, XJY, ix-1,iy  ,iz-1));
    mrc_fld_data_t cz = (F3S(p_f, XJZ, ix  ,iy  ,iz  ) +
		F3S(p_f, XJZ, ix-1,iy  ,iz  ) +
		F3S(p_f, XJZ, ix  ,iy-1,iz  ) +
		F3S(p_f, XJZ, ix-1,iy-1,iz  ));
    mrc_fld_data_t ffx = s2 * (cy * F3S(p_f, BZ, ix,iy,iz) -
		      cz * F3S(p_f, BY, ix,iy,iz));
    mrc_fld_data_t ffy = s2 * (cz * F3S(p_f, BX, ix,iy,iz) -
		      cx * F3S(p_f, BZ, ix,iy,iz));
    mrc_fld_data_t ffz = s2 * (cx * F3S(p_f, BY, ix,iy,iz) -
		      cy * F3S(p_f, BX, ix,iy,iz));
    mrc_fld_data_t duu = (ffx * F3S(p_f, _VX, ix,iy,iz) +
		 ffy * F3S(p_f, _VY, ix,iy,iz) +
		 ffz * F3S(p_f, _VZ, ix,iy,iz));

    F3S(p_f, m_next + _RV1X, ix,iy,iz) += ffx;
    F3S(p_f, m_next + _RV1Y, ix,iy,iz) += ffy;
    F3S(p_f, m_next + _RV1Z, ix,iy,iz) += ffz;
    F3S(p_f, m_next + _UU1 , ix,iy,iz) += duu;
  } fld3d_foreach_end;
}

static void
res1_const_c(fld3d_t p_f)
{
  mrc_fld_data_t diffsphere2 = sqr(s_diffsphere);

  fld3d_foreach(ix,iy,iz, 1, 1) {
    F3S(p_f, _RESIS, ix,iy,iz) = 0.f;
    mrc_fld_data_t r2 = FX2X(ix) + FX2Y(iy) + FX2Z(iz);
    if (r2 < diffsphere2)
      continue;
    if (iy + s_patch.off[1] < s_diff_obnd)
      continue;
    if (iz + s_patch.off[2] < s_diff_obnd)
      continue;
    if (ix + s_patch.off[0] >= s_gdims[0] - s_diff_obnd)
      continue;
    if (iy + s_patch.off[1] >= s_gdims[1] - s_diff_obnd)
      continue;
    if (iz + s_patch.off[2] >= s_gdims[2] - s_diff_obnd)
      continue;

    F3S(p_f, _RESIS, ix,iy,iz) = s_eta;
  } fld3d_foreach_end;
}

static void
calc_resis_const_c(fld3d_t p_f, int m_curr)
{
  curbc_c(p_f, m_curr);
  res1_const_c(p_f);
}

static void
calc_resis_nl1_c(fld3d_t p_f, int m_curr)
{
  // used to zero _RESIS field, but that's not needed.
}

static inline float
bcthy3f(mrc_fld_data_t s1, mrc_fld_data_t s2)
{
  if (s1 > 0.f && fabsf(s2) > REPS) {
/* .if(calce_aspect_low) then */
/* .call lowmask(IX, 0, 0,tl1) */
/* .call lowmask( 0,IY, 0,tl2) */
/* .call lowmask( 0, 0,IZ,tl3) */
/* .call lowmask(IX,IY,IZ,tl4) */
/*       tt=tt*(1.0-max(tl1,tl2,tl3,tl4)) */
    return s1 / s2;
  }
  return 0.f;
}

static inline void
calc_avg_dz_By(fld3d_t p_f, int m_curr, int XX, int YY, int ZZ,
	       int JX1, int JY1, int JZ1, int JX2, int JY2, int JZ2)
{
  // d_z B_y, d_y B_z on x edges
  fld3d_foreach(ix,iy,iz, 2, 1) {
    mrc_fld_data_t bd1[3] = { BD1X(ix), BD1Y(iy), BD1Z(iz) };

    F3S(p_f, _TMP1, ix,iy,iz) = bd1[ZZ] * 
      (F3S(p_f, m_curr + _B1X + YY, ix+JX2,iy+JY2,iz+JZ2) - F3S(p_f, m_curr + _B1X + YY, ix,iy,iz));
    F3S(p_f, _TMP2, ix,iy,iz) = bd1[YY] * 
      (F3S(p_f, m_curr + _B1X + ZZ, ix+JX1,iy+JY1,iz+JZ1) - F3S(p_f, m_curr + _B1X + ZZ, ix,iy,iz));
  } fld3d_foreach_end;

  // .5 * harmonic average if same sign
  fld3d_foreach(ix,iy,iz, 1, 1) {
    mrc_fld_data_t s1, s2;
    // dz_By on y face
    s1 = F3S(p_f, _TMP1, ix,iy,iz) * F3S(p_f, _TMP1, ix-JX2,iy-JY2,iz-JZ2);
    s2 = F3S(p_f, _TMP1, ix,iy,iz) + F3S(p_f, _TMP1, ix-JX2,iy-JY2,iz-JZ2);
    F3S(p_f, _TMP3, ix,iy,iz) = bcthy3f(s1, s2);
    // dy_Bz on z face
    s1 = F3S(p_f, _TMP2, ix,iy,iz) * F3S(p_f, _TMP2, ix-JX1,iy-JY1,iz-JZ1);
    s2 = F3S(p_f, _TMP2, ix,iy,iz) + F3S(p_f, _TMP2, ix-JX1,iy-JY1,iz-JZ1);
    F3S(p_f, _TMP4, ix,iy,iz) = bcthy3f(s1, s2);
  } fld3d_foreach_end;
}

#define CC_TO_EC(p_f, m, ix,iy,iz, IX,IY,IZ) \
  (.25f * (F3S(p_f, m, ix   ,iy   ,iz   ) +  \
	   F3S(p_f, m, ix   ,iy+IY,iz+IZ) +  \
	   F3S(p_f, m, ix+IX,iy   ,iz+IZ) +  \
	   F3S(p_f, m, ix+IX,iy+IY,iz   )))

static inline void
calc_v_x_B(mrc_fld_data_t ttmp[2], fld3d_t p_f, int m_curr, int ix, int iy, int iz,
	   int XX, int YY, int ZZ, int IX, int IY, int IZ,
	   int JX1, int JY1, int JZ1, int JX2, int JY2, int JZ2,
	   mrc_fld_data_t dt)
{
  mrc_fld_data_t bd2[3] = { BD2X(ix), BD2Y(iy), BD2Z(iz) };
  mrc_fld_data_t bd2p[3] = { BD2X(ix+1), BD2Y(iy+1), BD2Z(iz+1) };
  mrc_fld_data_t vbZZ;
  // edge centered velocity
  mrc_fld_data_t vvYY = CC_TO_EC(p_f, _VX + YY, ix,iy,iz, IX,IY,IZ) /* - d_i * vcurrYY */;
  if (vvYY > 0.f) {
    vbZZ = F3S(p_f, m_curr + _B1X + ZZ, ix,iy,iz) +
      F3S(p_f, _TMP4, ix,iy,iz) * (bd2[YY] - dt*vvYY);
  } else {
    vbZZ = F3S(p_f, m_curr + _B1X + ZZ, ix+JX1,iy+JY1,iz+JZ1) -
      F3S(p_f, _TMP4, ix+JX1,iy+JY1,iz+JZ1) * (bd2p[YY] + dt*vvYY);
  }
  ttmp[0] = vbZZ * vvYY;
  
  mrc_fld_data_t vbYY;
  // edge centered velocity
  mrc_fld_data_t vvZZ = CC_TO_EC(p_f, _VX + ZZ, ix,iy,iz, IX,IY,IZ) /* - d_i * vcurrZZ */;
  if (vvZZ > 0.f) {
    vbYY = F3S(p_f, m_curr + _B1X + YY, ix,iy,iz) +
      F3S(p_f, _TMP3, ix,iy,iz) * (bd2[ZZ] - dt*vvZZ);
  } else {
    vbYY = F3S(p_f, m_curr + _B1X + YY, ix+JX2,iy+JY2,iz+JZ2) -
      F3S(p_f, _TMP3, ix+JX2,iy+JY2,iz+JZ2) * (bd2p[ZZ] + dt*vvZZ);
  }
  ttmp[1] = vbYY * vvZZ;
}

static void
bcthy3z_NL1(fld3d_t p_f, int XX, int YY, int ZZ, int IX, int IY, int IZ,
	    int JX1, int JY1, int JZ1, int JX2, int JY2, int JZ2,
	    mrc_fld_data_t dt, int m_curr)
{
  calc_avg_dz_By(p_f, m_curr, XX, YY, ZZ, JX1, JY1, JZ1, JX2, JY2, JZ2);

  mrc_fld_data_t diffmul = 1.f;
  if (s_mhd_time < s_diff_timelo) { // no anomalous res at startup
    diffmul = 0.f;
  }

  // edge centered E = - v x B (+ dissipation)
  fld3d_foreach(ix,iy,iz, 1, 0) {
    mrc_fld_data_t ttmp[2];
    calc_v_x_B(ttmp, p_f, m_curr, ix, iy, iz, XX, YY, ZZ, IX, IY, IZ,
	       JX1, JY1, JZ1, JX2, JY2, JZ2, dt);
    
    mrc_fld_data_t t1m = F3S(p_f, m_curr + _B1X + ZZ, ix+JX1,iy+JY1,iz+JZ1) - F3S(p_f, m_curr + _B1X + ZZ, ix,iy,iz);
    mrc_fld_data_t t1p = fabsf(F3S(p_f, m_curr + _B1X + ZZ, ix+JX1,iy+JY1,iz+JZ1)) + fabsf(F3S(p_f, m_curr + _B1X + ZZ, ix,iy,iz));
    mrc_fld_data_t t2m = F3S(p_f, m_curr + _B1X + YY, ix+JX2,iy+JY2,iz+JZ2) - F3S(p_f, m_curr + _B1X + YY, ix,iy,iz);
    mrc_fld_data_t t2p = fabsf(F3S(p_f, m_curr + _B1X + YY, ix+JX2,iy+JY2,iz+JZ2)) + fabsf(F3S(p_f, m_curr + _B1X + YY, ix,iy,iz));
    mrc_fld_data_t tp = t1p + t2p + REPS;
    mrc_fld_data_t tpi = diffmul / tp;
    mrc_fld_data_t d1 = sqr(t1m * tpi);
    mrc_fld_data_t d2 = sqr(t2m * tpi);
    if (d1 < s_diffth) d1 = 0.;
    if (d2 < s_diffth) d2 = 0.;
    ttmp[0] -= d1 * t1m * F3S(p_f, _RMASK, ix,iy,iz);
    ttmp[1] -= d2 * t2m * F3S(p_f, _RMASK, ix,iy,iz);
    F3S(p_f, _RESIS, ix,iy,iz) += fabsf(d1+d2) * F3S(p_f, _ZMASK, ix,iy,iz);
    F3S(p_f, _FLX + XX, ix,iy,iz) = ttmp[0] - ttmp[1];
  } fld3d_foreach_end;
}

static void
bcthy3z_const(fld3d_t p_f, int XX, int YY, int ZZ, int IX, int IY, int IZ,
	      int JX1, int JY1, int JZ1, int JX2, int JY2, int JZ2, mrc_fld_data_t dt, int m_curr)
{
  calc_avg_dz_By(p_f, m_curr, XX, YY, ZZ, JX1, JY1, JZ1, JX2, JY2, JZ2);

  // edge centered E = - v x B (+ dissipation)
  fld3d_foreach(ix,iy,iz, 1, 0) {
    mrc_fld_data_t ttmp[2];
    calc_v_x_B(ttmp, p_f, m_curr, ix, iy, iz, XX, YY, ZZ, IX, IY, IZ,
	       JX1, JY1, JZ1, JX2, JY2, JZ2, dt);

    mrc_fld_data_t vcurrXX = CC_TO_EC(p_f, _CURRX + XX, ix,iy,iz, IX,IY,IZ);
    mrc_fld_data_t vresis = CC_TO_EC(p_f, _RESIS, ix,iy,iz, IX,IY,IZ);
    F3S(p_f, _FLX + XX, ix,iy,iz) = ttmp[0] - ttmp[1] - vresis * vcurrXX;
  } fld3d_foreach_end;
}

static void
calce_nl1_c(fld3d_t p_f, mrc_fld_data_t dt, int m_curr)
{
  bcthy3z_NL1(p_f, 0,1,2, 0,1,1, 0,1,0, 0,0,1, dt, m_curr);
  bcthy3z_NL1(p_f, 1,2,0, 1,0,1, 0,0,1, 1,0,0, dt, m_curr);
  bcthy3z_NL1(p_f, 2,0,1, 1,1,0, 1,0,0, 0,1,0, dt, m_curr);
}

static void
calce_const_c(fld3d_t p_f, mrc_fld_data_t dt, int m_curr)
{
  bcthy3z_const(p_f, 0,1,2, 0,1,1, 0,1,0, 0,0,1, dt, m_curr);
  bcthy3z_const(p_f, 1,2,0, 1,0,1, 0,0,1, 1,0,0, dt, m_curr);
  bcthy3z_const(p_f, 2,0,1, 1,1,0, 1,0,0, 0,1,0, dt, m_curr);
}

static void
calce_c(fld3d_t p_f, mrc_fld_data_t dt, int m_curr)
{
  switch (s_magdiffu) {
  case MAGDIFFU_NL1:
    return calce_nl1_c(p_f, dt, m_curr);
  case MAGDIFFU_CONST:
    return calce_const_c(p_f, dt, m_curr);
  default:
    assert(0);
  }
}

static void
bpush_c(fld3d_t p_f, mrc_fld_data_t dt, int m_prev, int m_next)
{
  fld3d_foreach(ix,iy,iz, 0, 0) {
    F3S(p_f, m_next + _B1X, ix,iy,iz) = F3S(p_f, m_prev + _B1X, ix,iy,iz) +
      dt * (BD3Y(iy) * (F3S(p_f,_FLZ, ix,iy,iz) - F3S(p_f,_FLZ, ix,iy-1,iz)) -
	    BD3Z(iz) * (F3S(p_f,_FLY, ix,iy,iz) - F3S(p_f,_FLY, ix,iy,iz-1)));
    F3S(p_f, m_next + _B1Y, ix,iy,iz) = F3S(p_f, m_prev + _B1Y, ix,iy,iz) +
      dt * (BD3Z(iz) * (F3S(p_f,_FLX, ix,iy,iz) - F3S(p_f,_FLX, ix,iy,iz-1)) -
	    BD3X(ix) * (F3S(p_f,_FLZ, ix,iy,iz) - F3S(p_f,_FLZ, ix-1,iy,iz)));
    F3S(p_f, m_next + _B1Z, ix,iy,iz) = F3S(p_f, m_prev + _B1Z, ix,iy,iz) +
      dt * (BD3X(ix) * (F3S(p_f,_FLY, ix,iy,iz) - F3S(p_f,_FLY, ix-1,iy,iz)) -
	    BD3Y(iy) * (F3S(p_f,_FLX, ix,iy,iz) - F3S(p_f,_FLX, ix,iy-1,iz)));
  } fld3d_foreach_end;
}

static void
pushstage_c(fld3d_t p_f, mrc_fld_data_t dt, int m_prev, int m_curr, int m_next,
	    int limit)
{
  rmaskn_c(p_f);

  if (limit != LIMIT_NONE) {
    vgrs(p_f, _BX, 0.f); vgrs(p_f, _BY, 0.f); vgrs(p_f, _BZ, 0.f);
    limit1_c(p_f, _PP, _BX);
    // limit2, 3
  }

  pushfv_c(p_f, _RR1 , dt, m_prev, m_curr, m_next, limit);
  pushfv_c(p_f, _RV1X, dt, m_prev, m_curr, m_next, limit);
  pushfv_c(p_f, _RV1Y, dt, m_prev, m_curr, m_next, limit);
  pushfv_c(p_f, _RV1Z, dt, m_prev, m_curr, m_next, limit);
  pushfv_c(p_f, _UU1 , dt, m_prev, m_curr, m_next, limit);

  pushpp_c(p_f, dt, m_next);

  switch (s_magdiffu) {
  case MAGDIFFU_NL1:
    calc_resis_nl1_c(p_f, m_curr);
    break;
  case MAGDIFFU_CONST:
    calc_resis_const_c(p_f, m_curr);
    break;
  default:
    assert(0);
  }

  push_ej_c(p_f, dt, m_curr, m_next);
  calce_c(p_f, dt, m_curr);
  bpush_c(p_f, dt, m_prev, m_next);
}

// ======================================================================
// ggcm_mhd_step subclass "c"
//
// this class will do full predictor / corrector steps,
// ie., including primvar() etc.

struct ggcm_mhd_step_c {
  struct mhd_options opt;
};

#define ggcm_mhd_step_c(step) mrc_to_subobj(step, struct ggcm_mhd_step_c)

// ----------------------------------------------------------------------
// ggcm_mhd_step_c_newstep

static void
ggcm_mhd_step_c_newstep(struct ggcm_mhd_step *step, float *p_dtn)
{
  struct ggcm_mhd *mhd = step->mhd;

  ggcm_mhd_fill_ghosts(mhd, mhd->fld, _RR1, mhd->time);

  fld3d_t p_f;
  fld3d_setup(&p_f, mhd->fld);

  static int PR;
  if (!PR) {
    PR = prof_register("newstep_sc_ggcm", 1., 0, 0);
  }
  prof_start(PR);

  mrc_fld_data_t splim2   = sqr(s_speedlimit_code);
  mrc_fld_data_t isphere2 = sqr(s_isphere);
  mrc_fld_data_t va02i    = 1.f / splim2;
  mrc_fld_data_t eps      = 1e-9f;
  mrc_fld_data_t epsz     = 1e-15f;
  mrc_fld_data_t gamma_m1 = s_gamma - 1.f;

  mrc_fld_data_t two_pi_d_i = 2. * M_PI * s_d_i;
  bool have_hall = s_d_i > 0.f;

  mrc_fld_data_t dt = 1e10f;
  pde_for_each_patch(p) {
    fld3d_get(&p_f, p);
    fld3d_foreach(i, j, k, 0, 0) {
      mrc_fld_data_t hh = mrc_fld_max(mrc_fld_max(FD1X(i), FD1Y(j)), FD1Z(k));
      mrc_fld_data_t rri = 1.f / mrc_fld_abs(F3S(p_f, RR, i,j,k)); // FIME abs necessary?
      mrc_fld_data_t bb = 
	sqr(.5f * (F3S(p_f, BX, i,j,k) + F3S(p_f, BX, i-1,j,k))) + 
	sqr(.5f * (F3S(p_f, BY, i,j,k) + F3S(p_f, BY, i,j-1,k))) +
	sqr(.5f * (F3S(p_f, BZ, i,j,k) + F3S(p_f, BZ, i,j,k-1)));
      if (have_hall) {
	bb *= 1 + sqr(two_pi_d_i * hh);
      }
      mrc_fld_data_t vv1 = fminf(bb * rri, splim2);
      
      mrc_fld_data_t rv2 = 
	sqr(F3S(p_f, RVX, i,j,k)) + sqr(F3S(p_f, RVY, i,j,k)) + sqr(F3S(p_f, RVZ, i,j,k));
      mrc_fld_data_t rvv = rri * rv2;
      mrc_fld_data_t pp = gamma_m1 * (F3S(p_f, UU, i,j,k) - .5f * rvv);
      mrc_fld_data_t vv2 = s_gamma * mrc_fld_max(0.f, pp) * rri;
      mrc_fld_data_t vv3 = rri * sqrtf(rv2);
      mrc_fld_data_t vv = sqrtf(vv1 + vv2) + vv3;
      vv = mrc_fld_max(eps, vv);
      
      mrc_fld_data_t ymask = 1.f;
      if (FX2X(i) + FX2Y(j) + FX2Z(k) < isphere2)
	ymask = 0.f;
      
      mrc_fld_data_t rrm = mrc_fld_max(epsz, bb * va02i);
      mrc_fld_data_t zmask = ymask * fminf(1.f, F3S(p_f, RR, i,j,k) / rrm);
      
      mrc_fld_data_t tt = s_cfl / mrc_fld_max(eps, hh*vv*zmask);
      dt = fminf(dt, tt);
    } fld3d_foreach_end;
    fld3d_put(&p_f, p);
  }
  mrc_fld_data_t dtn;
  MPI_Allreduce(&dt, &dtn, 1, MPI_MRC_FLD_DATA_T, MPI_MIN, ggcm_mhd_comm(mhd));
  
  mrc_fld_data_t dtmin = mhd->par.dtmin;
  // dtn is the global new timestep
  if (dtn <= dtmin) {
    mpi_printf(ggcm_mhd_comm(mhd), "!!! dt < dtmin. Dying now!\n");
    mpi_printf(ggcm_mhd_comm(mhd), "!!! dt %g -> %g, dtmin = %g\n",
	       mhd->dt, dtn, dtmin);   
    ggcm_mhd_wrongful_death(mhd, mhd->fld, -1);
  }

  prof_stop(PR);

  *p_dtn = dtn;
}

// ----------------------------------------------------------------------
// ggcm_mhd_step_c_pred

static void
ggcm_mhd_step_c_pred(struct ggcm_mhd_step *step)
{
  fld3d_t p_f;
  fld3d_setup(&p_f, step->mhd->fld);
  fld3d_get(&p_f, 0);
  pde_patch_set(0);
  s_mhd_time = step->mhd->time;

  patch_primvar(p_f, _RR1);
  patch_primbb(p_f, _RR1);
  patch_zmaskn(p_f);

  mrc_fld_data_t dth = .5f * step->mhd->dt;
  static int PR;
  if (!PR) {
    PR = prof_register("pred_c", 1., 0, 0);
  }
  prof_start(PR);
  pushstage_c(p_f, dth, _RR1, _RR1, _RR2, LIMIT_NONE);
  prof_stop(PR);

  fld3d_put(&p_f, 0);
}

// ----------------------------------------------------------------------
// ggcm_mhd_step_c_corr

static void
ggcm_mhd_step_c_corr(struct ggcm_mhd_step *step)
{
  fld3d_t p_f;
  fld3d_setup(&p_f, step->mhd->fld);
  fld3d_get(&p_f, 0);
  pde_patch_set(0);
  s_mhd_time = step->mhd->time;

  patch_primvar(p_f, _RR2);
  patch_primbb(p_f, _RR2);
  //  patch_zmaskn(p_f);

  static int PR;
  if (!PR) {
    PR = prof_register("corr_c", 1., 0, 0);
  }
  prof_start(PR);
  pushstage_c(p_f, step->mhd->dt, _RR1, _RR2, _RR1, LIMIT_1);
  prof_stop(PR);
  
  // --- check for NaNs and small density
  patch_badval_checks_sc(step->mhd, p_f, p_f);
  
  fld3d_put(&p_f, 0);
}

// ----------------------------------------------------------------------
// ggcm_mhd_step_c_setup_flds

static void
ggcm_mhd_step_c_setup_flds(struct ggcm_mhd_step *step)
{
  struct ggcm_mhd_step_c *sub = ggcm_mhd_step_c(step);
  struct ggcm_mhd *mhd = step->mhd;

  pde_mhd_set_options(mhd, &sub->opt);
  mrc_fld_set_type(mhd->fld, FLD_TYPE);
  mrc_fld_set_param_int(mhd->fld, "nr_ghosts", 2);
  mrc_fld_dict_add_int(mhd->fld, "mhd_type", MT_SEMI_CONSERVATIVE_GGCM);
  mrc_fld_set_param_int(mhd->fld, "nr_comps", _NR_FLDS);
}

// ----------------------------------------------------------------------
// ggcm_mhd_step_c_setup

static void
ggcm_mhd_step_c_setup(struct ggcm_mhd_step *step)
{
  struct ggcm_mhd *mhd = step->mhd;
  pde_setup(mhd->fld);
  pde_mhd_setup(mhd);

  mhd->ymask = mrc_fld_make_view(mhd->fld, _YMASK, _YMASK + 1);
  mrc_fld_set(mhd->ymask, 1.);

  ggcm_mhd_step_setup_member_objs_sub(step);
  ggcm_mhd_step_setup_super(step);

  s_fd1x = ggcm_mhd_crds_get_crd(mhd->crds, 0, FD1);
  s_fd1y = ggcm_mhd_crds_get_crd(mhd->crds, 1, FD1);
  s_fd1z = ggcm_mhd_crds_get_crd(mhd->crds, 2, FD1);
}

// ----------------------------------------------------------------------
// ggcm_mhd_step_c_destroy

static void
ggcm_mhd_step_c_destroy(struct ggcm_mhd_step *step)
{
  mrc_fld_destroy(step->mhd->ymask);
}

// ----------------------------------------------------------------------
// ggcm_mhd_step_c_get_e_ec

static void
ggcm_mhd_step_c_get_e_ec(struct ggcm_mhd_step *step, struct mrc_fld *Eout,
                         struct mrc_fld *state_vec)
{
  // the state vector should already be FLD_TYPE, but Eout is the data type
  // of the output
  struct mrc_fld *E = mrc_fld_get_as(Eout, FLD_TYPE);
  struct mrc_fld *x = mrc_fld_get_as(state_vec, FLD_TYPE);

  mrc_fld_foreach(x, ix, iy, iz, 1, 0) {
    F3(E, 0, ix, iy, iz) = F3(x, _FLX, ix, iy, iz);
    F3(E, 1, ix, iy, iz) = F3(x, _FLY, ix, iy, iz);
    F3(E, 2, ix, iy, iz) = F3(x, _FLZ, ix, iy, iz);
  } mrc_fld_foreach_end;

  mrc_fld_put_as(E, Eout);
  // FIXME, should use _put_as, but don't want copy-back
  if (strcmp(mrc_fld_type(state_vec), FLD_TYPE) != 0) {
    mrc_fld_destroy(x);
  }
} 

// ----------------------------------------------------------------------
// ggcm_mhd_step_c_diag_item_zmask_run

static void
ggcm_mhd_step_c_diag_item_zmask_run(struct ggcm_mhd_step *step,
				    struct ggcm_mhd_diag_item *item,
				    struct mrc_io *io, struct mrc_fld *f,
				    int diag_type, float plane)
{
  ggcm_mhd_diag_c_write_one_field(io, f, _ZMASK, "zmask", 1., diag_type, plane);
}

// ----------------------------------------------------------------------
// ggcm_mhd_step_c_diag_item_rmask_run

static void
ggcm_mhd_step_c_diag_item_rmask_run(struct ggcm_mhd_step *step,
				    struct ggcm_mhd_diag_item *item,
				    struct mrc_io *io, struct mrc_fld *f,
				    int diag_type, float plane)
{
  ggcm_mhd_diag_c_write_one_field(io, f, _RMASK, "rmask", 1., diag_type, plane);
}

// ----------------------------------------------------------------------
// subclass description

#define VAR(x) (void *)offsetof(struct ggcm_mhd_step_c, x)
static struct param ggcm_mhd_step_c_descr[] = {
  { "eqn"                , VAR(opt.eqn)            , PARAM_SELECT(OPT_EQN,
								  opt_eqn_descr)                },
  
  {},
};
#undef VAR

// ----------------------------------------------------------------------
// ggcm_mhd_step subclass "c_*"

struct ggcm_mhd_step_ops ggcm_mhd_step_c_ops = {
  .name                = ggcm_mhd_step_c_name,
  .size                = sizeof(struct ggcm_mhd_step_c),
  .param_descr         = ggcm_mhd_step_c_descr,
  .newstep             = ggcm_mhd_step_c_newstep,
  .pred                = ggcm_mhd_step_c_pred,
  .corr                = ggcm_mhd_step_c_corr,
  .run                 = ggcm_mhd_step_run_predcorr,
  .setup               = ggcm_mhd_step_c_setup,
  .destroy             = ggcm_mhd_step_c_destroy,
  .setup_flds          = ggcm_mhd_step_c_setup_flds,
  .get_e_ec            = ggcm_mhd_step_c_get_e_ec,
  .diag_item_zmask_run = ggcm_mhd_step_c_diag_item_zmask_run,
  .diag_item_rmask_run = ggcm_mhd_step_c_diag_item_rmask_run,
};
