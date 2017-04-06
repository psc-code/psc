
#ifndef PDE_MHD_PUSHFLUID_C
#define PDE_MHD_PUSHFLUID_C

static void
vgflrr_c(fld3d_t p_F, fld3d_t p_W)
{
  fld3d_foreach(i,j,k, 2, 2) {
    mrc_fld_data_t a = F3S(p_W, RR, i,j,k);
    F3S(p_F, 0, i,j,k) = a * F3S(p_W, VX, i,j,k);
    F3S(p_F, 1, i,j,k) = a * F3S(p_W, VY, i,j,k);
    F3S(p_F, 2, i,j,k) = a * F3S(p_W, VZ, i,j,k);
  } fld3d_foreach_end;
}

static void
vgflrvx_c(fld3d_t p_F, fld3d_t p_W)
{
  fld3d_foreach(i,j,k, 2, 2) {
    mrc_fld_data_t a = F3S(p_W, RR, i,j,k) * F3S(p_W, VX, i,j,k);
    F3S(p_F, 0, i,j,k) = a * F3S(p_W, VX, i,j,k);
    F3S(p_F, 1, i,j,k) = a * F3S(p_W, VY, i,j,k);
    F3S(p_F, 2, i,j,k) = a * F3S(p_W, VZ, i,j,k);
  } fld3d_foreach_end;
}

static void
vgflrvy_c(fld3d_t p_F, fld3d_t p_W)
{
  fld3d_foreach(i,j,k, 2, 2) {
    mrc_fld_data_t a = F3S(p_W, RR, i,j,k) * F3S(p_W, VY, i,j,k);
    F3S(p_F, 0, i,j,k) = a * F3S(p_W, VX, i,j,k);
    F3S(p_F, 1, i,j,k) = a * F3S(p_W, VY, i,j,k);
    F3S(p_F, 2, i,j,k) = a * F3S(p_W, VZ, i,j,k);
  } fld3d_foreach_end;
}

static void
vgflrvz_c(fld3d_t p_F, fld3d_t p_W)
{
  fld3d_foreach(i,j,k, 2, 2) {
    mrc_fld_data_t a = F3S(p_W, RR, i,j,k) * F3S(p_W, VZ, i,j,k);
    F3S(p_F, 0, i,j,k) = a * F3S(p_W, VX, i,j,k);
    F3S(p_F, 1, i,j,k) = a * F3S(p_W, VY, i,j,k);
    F3S(p_F, 2, i,j,k) = a * F3S(p_W, VZ, i,j,k);
  } fld3d_foreach_end;
}

static void
vgfluu_c(fld3d_t p_F, fld3d_t p_W)
{
  mrc_fld_data_t s = s_gamma / (s_gamma - 1.f);
  fld3d_foreach(i,j,k, 2, 2) {
    mrc_fld_data_t ep = s * F3S(p_W, PP, i,j,k) +
      .5f * F3S(p_W, RR, i,j,k) * (sqr(F3S(p_W, VX, i,j,k)) + 
				   sqr(F3S(p_W, VY, i,j,k)) + 
				   sqr(F3S(p_W, VZ, i,j,k)));
    F3S(p_F, 0, i,j,k) = ep * F3S(p_W, VX, i,j,k);
    F3S(p_F, 1, i,j,k) = ep * F3S(p_W, VY, i,j,k);
    F3S(p_F, 2, i,j,k) = ep * F3S(p_W, VZ, i,j,k);
  } fld3d_foreach_end;
}

static void
vgfl_c(fld3d_t p_F, fld3d_t p_W, int m)
{
  switch (m) {
  case RR:  return vgflrr_c(p_F, p_W);
  case RVX: return vgflrvx_c(p_F, p_W);
  case RVY: return vgflrvy_c(p_F, p_W);
  case RVZ: return vgflrvz_c(p_F, p_W);
  case UU:  return vgfluu_c(p_F, p_W);
  default: assert(0);
  }
}

static void
fluxl_c(fld3d_t p_Ffc, fld3d_t p_Fcc, fld3d_t p_cmsv, fld3d_t p_U, int m)
{
  fld3d_foreach(i,j,k, 1, 0) {
    mrc_fld_data_t aa = F3S(p_U, m, i,j,k);
    mrc_fld_data_t cmsv = F3S(p_cmsv, 0, i,j,k);
    F3S(p_Ffc, 0, i,j,k) =
      .5f * ((F3S(p_Fcc, 0, i,j,k) + F3S(p_Fcc, 0, i+di,j,k)) -
	     .5f * (F3S(p_cmsv, 0, i+di,j,k) + cmsv) * (F3S(p_U, m, i+di,j,k) - aa));
    F3S(p_Ffc, 1, i,j,k) =
      .5f * ((F3S(p_Fcc, 1, i,j,k) + F3S(p_Fcc, 1, i,j+dj,k)) -
	     .5f * (F3S(p_cmsv, 0, i,j+dj,k) + cmsv) * (F3S(p_U, m, i,j+dj,k) - aa));
    F3S(p_Ffc, 2, i,j,k) =
      .5f * ((F3S(p_Fcc, 2, i,j,k) + F3S(p_Fcc, 2, i,j,k+dk)) -
	     .5f * (F3S(p_cmsv, 0, i,j,k+dk) + cmsv) * (F3S(p_U, m, i,j,k+dk) - aa));
  } fld3d_foreach_end;
}

static void
fluxb_c(fld3d_t p_Ffc, fld3d_t p_Fcc, fld3d_t p_cmsv, fld3d_t p_U, int m, fld3d_t p_C)
{
  mrc_fld_data_t s1 = 1.f/12.f;
  mrc_fld_data_t s7 = 7.f * s1;

  fld3d_foreach(i,j,k, 1, 0) {
    mrc_fld_data_t fhx = (s7 * (F3S(p_Fcc, 0, i,j,k) + F3S(p_Fcc, 0, i+di,j,k)) -
		 s1 * (F3S(p_Fcc, 0, i-di,j,k) + F3S(p_Fcc, 0, i+2*di,j,k)));
    mrc_fld_data_t fhy = (s7 * (F3S(p_Fcc, 1, i,j,k) + F3S(p_Fcc, 1, i,j+dj,k)) -
		 s1 * (F3S(p_Fcc, 1, i,j-dj,k) + F3S(p_Fcc, 1, i,j+2*dj,k)));
    mrc_fld_data_t fhz = (s7 * (F3S(p_Fcc, 2, i,j,k) + F3S(p_Fcc, 2, i,j,k+dk)) -
		 s1 * (F3S(p_Fcc, 2, i,j,k-dk) + F3S(p_Fcc, 2, i,j,k+2*dk)));

    mrc_fld_data_t aa = F3S(p_U, m, i,j,k);
    mrc_fld_data_t cmsv = F3S(p_cmsv, 0, i,j,k);
    mrc_fld_data_t flx =
      .5f * ((F3S(p_Fcc, 0, i,j,k) + F3S(p_Fcc, 0, i+di,j,k)) -
	     .5f * (F3S(p_cmsv, 0, i+di,j,k) + cmsv) * (F3S(p_U, m, i+di,j,k) - aa));
    mrc_fld_data_t fly =
      .5f * ((F3S(p_Fcc, 1, i,j ,k) + F3S(p_Fcc, 1, i,j+dj,k)) -
	     .5f * (F3S(p_cmsv, 0, i,j+dj,k) + cmsv) * (F3S(p_U, m, i,j+dj,k) - aa));
    mrc_fld_data_t flz = 
      .5f * ((F3S(p_Fcc, 2, i,j,k) + F3S(p_Fcc, 2, i,j,k+dk)) -
	     .5f * (F3S(p_cmsv, 0, i,j,k+dk) + cmsv) * (F3S(p_U, m, i,j,k+dk) - aa));

    mrc_fld_data_t cx = F3S(p_C, 0, i,j,k);
    F3S(p_Ffc, 0, i,j,k) = cx * flx + (1.f - cx) * fhx;
    mrc_fld_data_t cy = F3S(p_C, 1, i,j,k);
    F3S(p_Ffc, 1, i,j,k) = cy * fly + (1.f - cy) * fhy;
    mrc_fld_data_t cz = F3S(p_C, 2, i,j,k);
    F3S(p_Ffc, 2, i,j,k) = cz * flz + (1.f - cz) * fhz;
  } fld3d_foreach_end;
}

static void
pushn_c(fld3d_t p_Unext, fld3d_t p_Uprev, fld3d_t p_F, fld3d_t p_ymask, int m, mrc_fld_data_t dt)
{
  if (p_Unext.arr_off == p_Uprev.arr_off) {
    fld3d_foreach(i,j,k, 0, 0) {
      mrc_fld_data_t s = dt * F3S(p_ymask, 0, i,j,k);
      F3S(p_Unext, m, i,j,k) += 
	- s * (FD1X(i) * (F3S(p_F, 0, i,j,k) - F3S(p_F, 0, i-di,j,k)) +
	       FD1Y(j) * (F3S(p_F, 1, i,j,k) - F3S(p_F, 1, i,j-dj,k)) +
	       FD1Z(k) * (F3S(p_F, 2, i,j,k) - F3S(p_F, 2, i,j,k-dk)));
    } fld3d_foreach_end;
  } else {
    fld3d_foreach(i,j,k, 0, 0) {
      mrc_fld_data_t s = dt * F3S(p_ymask, 0, i,j,k);
      F3S(p_Unext, m, i,j,k) = F3S(p_Uprev, m, i,j,k)
	- s * (FD1X(i) * (F3S(p_F, 0, i,j,k) - F3S(p_F, 0, i-di,j,k)) +
	       FD1Y(j) * (F3S(p_F, 1, i,j,k) - F3S(p_F, 1, i,j-dj,k)) +
	       FD1Z(k) * (F3S(p_F, 2, i,j,k) - F3S(p_F, 2, i,j,k-dk)));
    } fld3d_foreach_end;
  }
}

static void
vgrs(fld3d_t p_f, int m, mrc_fld_data_t s)
{
  fld3d_foreach(i,j,k, 2, 2) {
    F3S(p_f, m, i,j,k) = s;
  } fld3d_foreach_end;
}

static void
vgrv(fld3d_t p_to, int m_to, fld3d_t p_from, int m_from)
{
  fld3d_foreach(i,j,k, 2, 2) {
    F3S(p_to, m_to, i,j,k) = F3S(p_from, m_from, i,j,k);
  } fld3d_foreach_end;
}

static inline void
limit1a(fld3d_t p_U, int m, int i, int j, int k, int I, int J, int K, fld3d_t p_C, int C)
{
  const mrc_fld_data_t reps = 0.003;
  const mrc_fld_data_t seps = -0.001;
  const mrc_fld_data_t teps = 1.e-25;

  // Harten/Zwas type switch
  mrc_fld_data_t aa = F3S(p_U, m, i,j,k);
  mrc_fld_data_t a1 = F3S(p_U, m, i+I,j+J,k+K);
  mrc_fld_data_t a2 = F3S(p_U, m, i-I,j-J,k-K);
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
  F3S(p_C, C, i  ,j  ,k  ) = mrc_fld_max(F3S(p_C, C, i  ,j  ,k  ), r3);
  F3S(p_C, C, i-I,j-J,k-K) = mrc_fld_max(F3S(p_C, C, i-I,j-J,k-K), r3);
}

static void
limit1_c(fld3d_t p_U, int m, fld3d_t p_C)
{
  // we now don't limit at all if earlier than timelo
  /* if (s_mhd_time < s_timelo) { */
  /*   vgrs(p_C, 0, 1.f); */
  /*   vgrs(p_C, 1, 1.f); */
  /*   vgrs(p_C, 2, 1.f); */
  /*   return; */
  /* } */

  fld3d_foreach(i,j,k, 1, 1) {
    assert(!s_limit_aspect_low);
/* .if (limit_aspect_low) then */
/* .call lowmask(0,0,0,tl1) */
    limit1a(p_U, m, i,j,k, di,0,0, p_C, 0);
    limit1a(p_U, m, i,j,k, 0,dj,0, p_C, 1);
    limit1a(p_U, m, i,j,k, 0,0,dk, p_C, 2);
  } fld3d_foreach_end;
}

static void
pushfv_c(fld3d_t p_Unext, fld3d_t p_Uprev, fld3d_t p_Ucurr, int m, 
	 fld3d_t p_Wcurr, fld3d_t p_cmsv, fld3d_t p_ymask, mrc_fld_data_t dt,
	 bool limit, fld3d_t p_B)
{
  static fld3d_t p_Ffc, p_Fcc, p_C;
  fld3d_setup_tmp_compat(&p_Ffc, 3, _FLX);
  fld3d_setup_tmp_compat(&p_Fcc, 3, _TMP1);
  fld3d_setup_tmp_compat(&p_C, 3, _CX);

  vgfl_c(p_Fcc, p_Wcurr, m);
  if (!limit) {
    fluxl_c(p_Ffc, p_Fcc, p_cmsv, p_Ucurr, m);
  } else {
    vgrv(p_C, 0, p_B, 0);
    vgrv(p_C, 1, p_B, 1);
    vgrv(p_C, 2, p_B, 2);
    limit1_c(p_Ucurr, m, p_C);
    fluxb_c(p_Ffc, p_Fcc, p_cmsv, p_Ucurr, m, p_C);
  }

  pushn_c(p_Unext, p_Uprev, p_Ffc, p_ymask, m, dt);
}

static void
pushpp_c(fld3d_t p_Unext, fld3d_t p_W, fld3d_t p_zmask, mrc_fld_data_t dt)
{
  mrc_fld_data_t dth = -.5f * dt;
  fld3d_foreach(i,j,k, 0, 0) {
    mrc_fld_data_t z = dth * F3S(p_zmask, 0, i,j,k);
    F3S(p_Unext, RVX, i,j,k) += z * FD1X(i) * (F3S(p_W, PP, i+di,j,k) - F3S(p_W, PP, i-di,j,k));
    F3S(p_Unext, RVY, i,j,k) += z * FD1Y(j) * (F3S(p_W, PP, i,j+dj,k) - F3S(p_W, PP, i,j-dj,k));
    F3S(p_Unext, RVZ, i,j,k) += z * FD1Z(k) * (F3S(p_W, PP, i,j,k+dk) - F3S(p_W, PP, i,j,k-dk));
  } fld3d_foreach_end;
}

static void
patch_pushfluid_c(fld3d_t p_Unext, mrc_fld_data_t dt, fld3d_t p_Uprev,
		  fld3d_t p_Ucurr, fld3d_t p_W, fld3d_t p_cmsv, fld3d_t p_ymask,
		  fld3d_t p_zmask, int stage)
{
  static fld3d_t p_B;
  fld3d_setup_tmp_compat(&p_B, 3, _BX);
  bool limit = stage != 0 && s_mhd_time > s_timelo;

  if (limit) {
    vgrs(p_B, 0, 0.f); vgrs(p_B, 1, 0.f); vgrs(p_B, 2, 0.f);
    assert(!s_do_limit2);
    assert(!s_do_limit3);
    limit1_c(p_W, PP, p_B);
  }

  pushfv_c(p_Unext, p_Uprev, p_Ucurr, RR , p_W, p_cmsv, p_ymask, dt, limit, p_B);
  pushfv_c(p_Unext, p_Uprev, p_Ucurr, RVX, p_W, p_cmsv, p_ymask, dt, limit, p_B);
  pushfv_c(p_Unext, p_Uprev, p_Ucurr, RVY, p_W, p_cmsv, p_ymask, dt, limit, p_B);
  pushfv_c(p_Unext, p_Uprev, p_Ucurr, RVZ, p_W, p_cmsv, p_ymask, dt, limit, p_B);
  pushfv_c(p_Unext, p_Uprev, p_Ucurr, UU , p_W, p_cmsv, p_ymask, dt, limit, p_B);

  pushpp_c(p_Unext, p_W, p_zmask, dt);
}

// ----------------------------------------------------------------------
// patch_pushfluid_fortran

#if defined(HAVE_OPENGGCM_FORTRAN) && defined(MRC_FLD_AS_FLOAT_H)

#include "pde/pde_fortran.h"

#define pushfluid1_F77 F77_FUNC(pushfluid1,PUSHFLUID1)
#define pushfluid2_F77 F77_FUNC(pushfluid2,PUSHFLUID2)


void pushfluid1_F77(real *rr1, real *rv1x, real *rv1y, real *rv1z, real *uu1,
		    real *rr2, real *rv2x, real *rv2y, real *rv2z, real *uu2,
		    real *rr, real *vx, real *vy, real *vz, real *pp,
		    real *ymask, real *zmask, real *cmsv,
		    real *dth);

void pushfluid2_F77(real *rr1, real *rv1x, real *rv1y, real *rv1z, real *uu1,
		    real *rr2, real *rv2x, real *rv2y, real *rv2z, real *uu2,
		    real *rr, real *vx, real *vy, real *vz, real *pp,
		    real *ymask, real *zmask, real *cmsv,
		    real *dth, real *time);

static void
patch_pushfluid1_fortran(mrc_fld_data_t dth)
{
  pushfluid1_F77(F(s_p_f, _RR1), F(s_p_f, _RV1X), F(s_p_f, _RV1Y), F(s_p_f, _RV1Z), F(s_p_f, _UU1),
		 F(s_p_f, _RR2), F(s_p_f, _RV2X), F(s_p_f, _RV2Y), F(s_p_f, _RV2Z), F(s_p_f, _UU2),
		 F(s_p_f, _RR), F(s_p_f, _VX), F(s_p_f, _VY), F(s_p_f, _VZ), F(s_p_f, _PP),
		 F(s_p_f, _YMASK), F(s_p_f, _ZMASK), F(s_p_f, _CMSV),
		 &dth);
}

static void
patch_pushfluid2_fortran(mrc_fld_data_t dth)
{
  pushfluid2_F77(F(s_p_f, _RR1), F(s_p_f, _RV1X), F(s_p_f, _RV1Y), F(s_p_f, _RV1Z), F(s_p_f, _UU1),
		 F(s_p_f, _RR2), F(s_p_f, _RV2X), F(s_p_f, _RV2Y), F(s_p_f, _RV2Z), F(s_p_f, _UU2),
		 F(s_p_f, _RR), F(s_p_f, _VX), F(s_p_f, _VY), F(s_p_f, _VZ), F(s_p_f, _PP),
		 F(s_p_f, _YMASK), F(s_p_f, _ZMASK), F(s_p_f, _CMSV),
		 &dth, &s_mhd_time);
}

static void
patch_pushfluid_fortran(mrc_fld_data_t dt, int stage)
{
  if (stage == 0) {
    patch_pushfluid1_fortran(dt);
  } else {
    patch_pushfluid2_fortran(dt);
  }
}

#endif

// ----------------------------------------------------------------------
// patch_pushfluid

static void _mrc_unused
patch_pushfluid(fld3d_t p_Unext, mrc_fld_data_t dt, fld3d_t p_Uprev,
		fld3d_t p_Ucurr, fld3d_t p_W, fld3d_t p_cmsv, fld3d_t p_ymask,
		fld3d_t p_zmask, int stage)
{
  int opt_mhd_pushfluid = stage ? s_opt_mhd_pushfluid2 : s_opt_mhd_pushfluid1;

  if (opt_mhd_pushfluid == OPT_MHD_C) {
    patch_pushfluid_c(p_Unext, dt, p_Uprev, p_Ucurr, p_W,
		      p_cmsv, p_ymask, p_zmask, stage);
#if defined(HAVE_OPENGGCM_FORTRAN) && defined(MRC_FLD_AS_FLOAT_H)
  } else if (opt_mhd_pushfluid == OPT_MHD_FORTRAN) {
    patch_pushfluid_fortran(dt, stage);
#endif
  } else {
    assert(0);
  }
}

#endif
