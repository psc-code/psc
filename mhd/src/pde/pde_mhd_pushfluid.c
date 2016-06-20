
#ifndef PDE_MHD_PUSHFLUID_C
#define PDE_MHD_PUSHFLUID_C

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
    assert(!s_limit_aspect_low);
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
pushfluid_c(fld3d_t p_f, mrc_fld_data_t dt,
	    int m_prev, int m_curr, int m_next, int limit)
{
  if (limit != LIMIT_NONE) {
    vgrs(p_f, _BX, 0.f); vgrs(p_f, _BY, 0.f); vgrs(p_f, _BZ, 0.f);
    assert(!s_do_limit2);
    assert(!s_do_limit3);
    limit1_c(p_f, _PP, _BX);
  }

  pushfv_c(p_f, _RR1 , dt, m_prev, m_curr, m_next, limit);
  pushfv_c(p_f, _RV1X, dt, m_prev, m_curr, m_next, limit);
  pushfv_c(p_f, _RV1Y, dt, m_prev, m_curr, m_next, limit);
  pushfv_c(p_f, _RV1Z, dt, m_prev, m_curr, m_next, limit);
  pushfv_c(p_f, _UU1 , dt, m_prev, m_curr, m_next, limit);

  pushpp_c(p_f, dt, m_next);
}

// ----------------------------------------------------------------------
// patch_pushfluid1_c

static void
patch_pushfluid1_c(fld3d_t p_f, mrc_fld_data_t dt)
{
  pushfluid_c(p_f, dt, _RR1, _RR1, _RR2, LIMIT_NONE);
}

// ----------------------------------------------------------------------
// patch_pushfluid1_fortran

#if defined(HAVE_OPENGGCM_FORTRAN) && defined(MRC_FLD_AS_FLOAT_H)

#include "pde/pde_fortran.h"

#define pushfluid1_F77 F77_FUNC(pushfluid1,PUSHFLUID1)

void pushfluid1_F77(real *rr1, real *rv1x, real *rv1y, real *rv1z, real *uu1,
		    real *rr2, real *rv2x, real *rv2y, real *rv2z, real *uu2,
		    real *rr, real *vx, real *vy, real *vz, real *pp,
		    real *ymask, real *zmask, real *cmsv,
		    real *dth);

static void
patch_pushfluid1_fortran(fld3d_t p_f, mrc_fld_data_t dth)
{
  pushfluid1_F77(F(p_f, _RR1), F(p_f, _RV1X), F(p_f, _RV1Y), F(p_f, _RV1Z), F(p_f, _UU1),
		 F(p_f, _RR2), F(p_f, _RV2X), F(p_f, _RV2Y), F(p_f, _RV2Z), F(p_f, _UU2),
		 F(p_f, _RR), F(p_f, _VX), F(p_f, _VY), F(p_f, _VZ), F(p_f, _PP),
		 F(p_f, _YMASK), F(p_f, _ZMASK), F(p_f, _CMSV),
		 &dth);
}

#endif

// ----------------------------------------------------------------------
// patch_pushfluid1

static void _mrc_unused
patch_pushfluid1(fld3d_t p_f, mrc_fld_data_t dt)
{
  if (s_opt_mhd_pushfluid1 == OPT_MHD_C) {
    patch_pushfluid1_c(p_f, dt);
#if defined(HAVE_OPENGGCM_FORTRAN) && defined(MRC_FLD_AS_FLOAT_H)
  } else if (s_opt_mhd_pushfluid1 == OPT_MHD_FORTRAN) {
    patch_pushfluid1_fortran(p_f, dt);
#endif
  } else {
    assert(0);
  }
}

// ======================================================================

// ----------------------------------------------------------------------
// patch_pushfluid2_c

static void
patch_pushfluid2_c(fld3d_t p_f, mrc_fld_data_t dt)
{
  pushfluid_c(p_f, dt, _RR1, _RR2, _RR1, LIMIT_1);
}

// ----------------------------------------------------------------------
// patch_pushfluid2_fortran

#if defined(HAVE_OPENGGCM_FORTRAN) && defined(MRC_FLD_AS_FLOAT_H)

#include "pde/pde_fortran.h"

#define pushfluid2_F77 F77_FUNC(pushfluid2,PUSHFLUID2)

void pushfluid2_F77(real *rr1, real *rv1x, real *rv1y, real *rv1z, real *uu1,
		    real *rr2, real *rv2x, real *rv2y, real *rv2z, real *uu2,
		    real *rr, real *vx, real *vy, real *vz, real *pp,
		    real *ymask, real *zmask, real *cmsv,
		    real *dth, real *time);

static void
patch_pushfluid2_fortran(fld3d_t p_f, mrc_fld_data_t dth)
{
  pushfluid2_F77(F(p_f, _RR1), F(p_f, _RV1X), F(p_f, _RV1Y), F(p_f, _RV1Z), F(p_f, _UU1),
		 F(p_f, _RR2), F(p_f, _RV2X), F(p_f, _RV2Y), F(p_f, _RV2Z), F(p_f, _UU2),
		 F(p_f, _RR), F(p_f, _VX), F(p_f, _VY), F(p_f, _VZ), F(p_f, _PP),
		 F(p_f, _YMASK), F(p_f, _ZMASK), F(p_f, _CMSV),
		 &dth, &s_mhd_time);
}

#endif

// ----------------------------------------------------------------------
// patch_pushfluid2

static void _mrc_unused
patch_pushfluid2(fld3d_t p_f, mrc_fld_data_t dt)
{
  if (s_opt_mhd_pushfluid2 == OPT_MHD_C) {
    patch_pushfluid2_c(p_f, dt);
#if defined(HAVE_OPENGGCM_FORTRAN) && defined(MRC_FLD_AS_FLOAT_H)
  } else if (s_opt_mhd_pushfluid2 == OPT_MHD_FORTRAN) {
    patch_pushfluid2_fortran(p_f, dt);
#endif
  } else {
    assert(0);
  }
}

#endif
