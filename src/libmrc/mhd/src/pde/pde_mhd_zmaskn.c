
#ifndef PDE_MHD_ZMASKN_C
#define PDE_MHD_ZMASKN_C

#include "pde/pde_mhd_primbb.c"

// ----------------------------------------------------------------------
// patch_calc_zmask_c

static void
patch_calc_zmask_c(fld3d_t p_zmask, fld3d_t p_U, fld3d_t p_ymask)
{
  fld3d_t p_B = fld3d_make_view(p_U, BX);

  mrc_fld_data_t va02i = 1.f / sqr(s_speedlimit_code);
  mrc_fld_data_t eps = 1e-30f;

  fld3d_foreach(i,j,k, 1, 1) {
    float bb = (sqr(BTXcc(p_B, i,j,k)) + sqr(BTYcc(p_B, i,j,k)) + sqr(BTZcc(p_B, i,j,k)));
    float rrm = mrc_fld_max(eps, s_mu0_inv * bb * va02i);
    F3S(p_zmask, 0, i,j,k) = F3S(p_ymask, 0, i,j,k) * 
      mrc_fld_min(1.f, F3S(p_U, RR, i,j,k) / rrm);
  } fld3d_foreach_end;
}

// ----------------------------------------------------------------------
// patch_zmaskn_c
//
// like the above, but expects pre-calculated bcc
// (superseded)

static void
patch_zmaskn_c(fld3d_t p_zmask, fld3d_t p_W, fld3d_t p_bcc, fld3d_t p_ymask)
{
  mrc_fld_data_t va02i = 1.f / sqr(s_speedlimit_code);
  mrc_fld_data_t eps = 1e-30f;

  fld3d_foreach(i,j,k, 2, 2) {
    float bb = (sqr(F3S(p_bcc, 0, i,j,k)) + 
		sqr(F3S(p_bcc, 1, i,j,k)) +
		sqr(F3S(p_bcc, 2, i,j,k)));
    float rrm = mrc_fld_max(eps, s_mu0_inv * bb * va02i);
    F3S(p_zmask, 0, i,j,k) = F3S(p_ymask, 0, i,j,k) * 
      mrc_fld_min(1.f, F3S(p_W, RR, i,j,k) / rrm);
  } fld3d_foreach_end;
}

// ----------------------------------------------------------------------
// patch_zmaskn_fortran

#if defined(HAVE_OPENGGCM_FORTRAN) && defined(MRC_FLD_AS_FLOAT_H)

#include "pde/pde_fortran.h"

#define zmaskn_F77 F77_FUNC(zmaskn,ZMASKN)

void zmaskn_F77(real *rr, real *bx, real *by, real *bz, 
		real *zmask, real *ymask);

static void
patch_zmaskn_fortran(fld3d_t p_zmask, fld3d_t p_W, fld3d_t p_bcc, fld3d_t p_ymask)
{
  zmaskn_F77(F(p_W, RR), F(p_bcc, 0), F(p_bcc, 1), F(p_bcc, 2),
	     F(p_zmask, 0), F(p_ymask, 0));
}

#endif

// ----------------------------------------------------------------------
// patch_zmaskn

static void _mrc_unused
patch_zmaskn(fld3d_t p_zmask, fld3d_t p_W, fld3d_t p_bcc, fld3d_t p_ymask)
{
  if (s_opt_mhd_zmaskn == OPT_MHD_C) {
    patch_zmaskn_c(p_zmask, p_W, p_bcc, p_ymask);
#if defined(HAVE_OPENGGCM_FORTRAN) && defined(MRC_FLD_AS_FLOAT_H)
  } else if (s_opt_mhd_zmaskn == OPT_MHD_FORTRAN) {
    patch_zmaskn_fortran(p_zmask, p_W, p_bcc, p_ymask);
#endif
  } else {
    assert(0);
  }
}

// ----------------------------------------------------------------------
// patch_calc_zmask_gold

static void _mrc_unused
patch_calc_zmask_gold(fld3d_t p_zmask, fld3d_t p_U, fld3d_t p_ymask)
{
  static fld3d_t p_bcc;
  fld3d_setup_tmp_compat(&p_bcc, 3, _BX);

  patch_primbb(p_bcc, p_U);
  // we kinda should be passing p_W rather than p_U to zmaskn(),
  // but since we only access RR, either works
  patch_zmaskn(p_zmask, p_U, p_bcc, p_ymask);
}

// ----------------------------------------------------------------------
// patch_calc_zmask
//
// like zmaskn(), but doesn't require preparation step (primbb)

static void _mrc_unused
patch_calc_zmask(fld3d_t p_zmask, fld3d_t p_U, fld3d_t p_ymask)
{
  patch_calc_zmask_c(p_zmask, p_U, p_ymask);
}

#endif
