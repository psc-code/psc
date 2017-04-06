
#ifndef PDE_MHD_PRIMVAR_C
#define PDE_MHD_PRIMVAR_C

#include "pde/pde_mhd_convert.c"

// ----------------------------------------------------------------------
// patch_cmsv

static void
patch_cmsv(fld3d_t p_cmsv, fld3d_t p_W, fld3d_t p_U)
{
  fld3d_foreach(i,j,k, 2, 2) {
    mrc_fld_data_t rri = 1.f / F3S(p_U, RR, i,j,k);
    mrc_fld_data_t rvv =
      F3S(p_W, VX, i,j,k) * F3S(p_U, RVX, i,j,k) +
      F3S(p_W, VY, i,j,k) * F3S(p_U, RVY, i,j,k) +
      F3S(p_W, VZ, i,j,k) * F3S(p_U, RVZ, i,j,k);
    mrc_fld_data_t cs2 = mrc_fld_max(s_gamma * F3S(p_W, PP, i,j,k) * rri, 0.f);
    F3S(p_cmsv, 0, i,j,k) = mrc_fld_sqrt(rvv * rri) + mrc_fld_sqrt(cs2);
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
patch_primvar_fortran(fld3d_t p_W, fld3d_t p_U, fld3d_t p_cmsv)
{
  int ntot = s_lgdims[0] * s_lgdims[1] * s_lgdims[2];
  // very ugly way to pass _TMP1...
  primvar_F77(&ntot, F(p_U, RR), F(p_U, RVX), F(p_U, RVY), F(p_U, RVZ), F(p_U, UU),
	      F(p_W, RR), F(p_W, VX), F(p_W, VY), F(p_W, VZ), F(p_W, PP), F(p_cmsv, 0),
	      &s_gamma, F(p_W, _TMP1 - _RR));
}

#endif

// ----------------------------------------------------------------------
// patch_primvar

static void _mrc_unused
patch_primvar(fld3d_t p_W, fld3d_t p_U, fld3d_t p_cmsv)
{
  if (s_opt_mhd_primvar == OPT_MHD_C) {
    patch_prim_from_cons_v2(p_W, p_U, 2);
    patch_cmsv(p_cmsv, p_W, p_U);
#if defined(HAVE_OPENGGCM_FORTRAN) && defined(MRC_FLD_AS_FLOAT_H)
  } else if (s_opt_mhd_primvar == OPT_MHD_FORTRAN) {
    patch_primvar_fortran(p_W, p_U, p_cmsv);
#endif
  } else {
    assert(0);
  }
}

#endif
