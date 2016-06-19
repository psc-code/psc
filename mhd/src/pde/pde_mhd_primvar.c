
// ----------------------------------------------------------------------
// patch_prim_from_cons_v2

static void
patch_prim_from_cons_v2(fld3d_t p_W, fld3d_t p_U, int sw)
{
  mrc_fld_data_t gamma_m1 = s_gamma - 1.f;

  fld3d_foreach(i,j,k, sw, sw) {
    F3S(p_W, RR, i,j,k) = F3S(p_U, RR, i,j,k);
    mrc_fld_data_t rri  = 1.f / F3S(p_U, RR, i,j,k);
    F3S(p_W, VX, i,j,k) = rri * F3S(p_U, RVX, i,j,k);
    F3S(p_W, VY, i,j,k) = rri * F3S(p_U, RVY, i,j,k);
    F3S(p_W, VZ, i,j,k) = rri * F3S(p_U, RVZ, i,j,k);
    mrc_fld_data_t rvv =
      F3S(p_W, VX, i,j,k) * F3S(p_U, RVX, i,j,k) +
      F3S(p_W, VY, i,j,k) * F3S(p_U, RVY, i,j,k) +
      F3S(p_W, VZ, i,j,k) * F3S(p_U, RVZ, i,j,k);
    F3S(p_W, PP, i,j,k) = gamma_m1 * (F3S(p_U, UU, i,j,k) - .5f * rvv);
  } fld3d_foreach_end;
}

// ----------------------------------------------------------------------
// patch_primvar_c

static void
patch_primvar_c(fld3d_t p_f, int m)
{
  fld3d_t p_W, p_U;
  fld3d_setup_view(&p_W, p_f, _RR);
  fld3d_setup_view(&p_U, p_f, m);
  patch_prim_from_cons_v2(p_W, p_U, 2);

  fld3d_foreach(i,j,k, 2, 2) {
    mrc_fld_data_t rri = 1.f / F3S(p_f, m + _RR1, i,j,k);
    mrc_fld_data_t rvv =
      F3S(p_f,_VX, i,j,k) * F3S(p_f, m + _RV1X, i,j,k) +
      F3S(p_f,_VY, i,j,k) * F3S(p_f, m + _RV1Y, i,j,k) +
      F3S(p_f,_VZ, i,j,k) * F3S(p_f, m + _RV1Z, i,j,k);
    mrc_fld_data_t cs2 = mrc_fld_max(s_gamma * F3S(p_f,_PP, i,j,k) * rri, 0.f);
    F3S(p_f,_CMSV, i,j,k) = mrc_fld_sqrt(rvv * rri) + mrc_fld_sqrt(cs2);
  } fld3d_foreach_end;
}

// ----------------------------------------------------------------------
// patch_primvar_fortran

#if defined(HAVE_OPENGGCM_FORTRAN) && defined(MRC_FLD_AS_FLOAT_H)

#include "pde/pde_fortran.h"

#define primvar_F77 F77_FUNC(primvar,PRIMVAR)

void primvar_F77(integer *ntot, 
		 real *rr1, real *rvx1, real *rvy1, real *rvz1, 
		 real *uu1, real *rr, real  *vx, real *vy, real *vz, real *pp,
		 real *cmsv, real *gamma, real *tmp1);

static void
patch_primvar_fortran(fld3d_t p_f, int m)
{
  int ntot = s_lgdims[0] * s_lgdims[1] * s_lgdims[2];
  if (m == _RR1) {
    primvar_F77(&ntot, F(p_f, _RR1), F(p_f, _RV1X), F(p_f, _RV1Y), F(p_f, _RV1Z), F(p_f, _UU1),
		F(p_f, _RR), F(p_f, _VX), F(p_f, _VY), F(p_f, _VZ), F(p_f, _PP), F(p_f, _CMSV),
		&s_gamma, F(p_f, _TMP1));
  } else if (m == _RR2) {
    primvar_F77(&ntot, F(p_f, _RR2), F(p_f, _RV2X), F(p_f, _RV2Y), F(p_f, _RV2Z), F(p_f, _UU2),
		F(p_f, _RR), F(p_f, _VX), F(p_f, _VY), F(p_f, _VZ), F(p_f, _PP), F(p_f, _CMSV),
		&s_gamma, F(p_f, _TMP1));
  } else {
    assert(0);
  }
}

#endif

// ----------------------------------------------------------------------
// patch_primvar

static void _mrc_unused
patch_primvar(fld3d_t p_f, int m)
{
  if (s_opt_mhd_primvar == OPT_MHD_C) {
    patch_primvar_c(p_f, m);
#if defined(HAVE_OPENGGCM_FORTRAN) && defined(MRC_FLD_AS_FLOAT_H)
  } else if (s_opt_mhd_primvar == OPT_MHD_FORTRAN) {
    patch_primvar_fortran(p_f, m);
#endif
  } else {
    assert(0);
  }
}

