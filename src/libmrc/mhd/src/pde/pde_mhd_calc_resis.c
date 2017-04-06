
#ifndef PDE_MHD_CALC_RESIS_C
#define PDE_MHD_CALC_RESIS_C

#include "pde/pde_mhd_calc_current.c"

// ======================================================================
// nl1

// ----------------------------------------------------------------------
// patch_calc_resis_nl1_c

static void
patch_calc_resis_nl1_c(fld3d_t p_resis)
{
  // FIXME, unnecessary
  fld3d_foreach(i,j,k, 2, 2) {
    F3S(p_resis, 0, i,j,k) = 0.f;
  } fld3d_foreach_end;
}

// ======================================================================
// const

// ----------------------------------------------------------------------
// patch_res1_const

static void _mrc_unused
patch_res1_const(fld3d_t p_resis)
{
  mrc_fld_data_t diffsphere2 = sqr(s_diffsphere);

  fld3d_foreach(i,j,k, 1, 1) {
    F3S(p_resis, 0, i,j,k) = 0.f;
    mrc_fld_data_t r2 = sqr(PDE_CRDX_CC(i)) + sqr(PDE_CRDY_CC(j)) + sqr(PDE_CRDZ_CC(k));
    if (r2 < diffsphere2)
      continue;
    if (j + s_patch.off[1] < s_diff_obnd)
      continue;
    if (k + s_patch.off[2] < s_diff_obnd)
      continue;
    if (i + s_patch.off[0] >= s_gdims[0] - s_diff_obnd)
      continue;
    if (j + s_patch.off[1] >= s_gdims[1] - s_diff_obnd)
      continue;
    if (k + s_patch.off[2] >= s_gdims[2] - s_diff_obnd)
      continue;

    F3S(p_resis, 0, i,j,k) = s_eta;
  } fld3d_foreach_end;
}

// ----------------------------------------------------------------------
// patch_calc_resis_const_c

static void
patch_calc_resis_const_c(fld3d_t p_resis, fld3d_t p_Jcc, fld3d_t p_U,
			 fld3d_t p_zmask)
{
  patch_calc_current_cc(p_Jcc, p_U, p_zmask);
  patch_res1_const(p_resis);
}

// ======================================================================
// fortran

// ----------------------------------------------------------------------
// patch_calc_resis_nl1/const_fortran

#if defined(HAVE_OPENGGCM_FORTRAN) && defined(MRC_FLD_AS_FLOAT_H)

#define calc_resis_nl1_F77 F77_FUNC(calc_resis_nl1,CALC_RESIS_NL1)
#define calc_resis_const_F77 F77_FUNC(calc_resis_const,CALC_RESIS_CONST)

void calc_resis_nl1_F77(real *bx, real *by, real *bz, real *resis);
void calc_resis_const_F77(real *bx, real *by, real *bz, 
			  real *currx, real *curry, real *currz,
			  real *tmp1, real *tmp2, real *tmp3,
			  real *flx, real *fly, real *flz,
			  real *zmask, real *rr, real *pp, real *resis);

static void
patch_calc_resis_nl1_fortran(fld3d_t p_resis)
{
  calc_resis_nl1_F77(NULL, NULL, NULL, F(p_resis, 0));
}

static void
patch_calc_resis_const_fortran(fld3d_t p_resis, fld3d_t p_Jcc, fld3d_t p_U, 
			       fld3d_t p_zmask)
{
  calc_resis_const_F77(F(p_U, BX), F(p_U, BY), F(p_U, BZ),
		       F(p_Jcc, 0), F(p_Jcc, 1), F(p_Jcc, 2),
		       F(s_p_f, _TMP1), F(s_p_f, _TMP2), F(s_p_f, _TMP3),
		       F(s_p_f, _FLX), F(s_p_f, _FLY), F(s_p_f, _FLZ),
		       F(p_zmask, 0), NULL, NULL, F(p_resis, 0));
}

#endif

// ======================================================================
// dispatch

// ----------------------------------------------------------------------
// patch_calc_resis_nl1

static void _mrc_unused
patch_calc_resis_nl1(fld3d_t p_resis)
{
  if (s_opt_mhd_calc_resis == OPT_MHD_C) {
    patch_calc_resis_nl1_c(p_resis);
#if defined(HAVE_OPENGGCM_FORTRAN) && defined(MRC_FLD_AS_FLOAT_H)
  } else if (s_opt_mhd_calc_resis == OPT_MHD_FORTRAN) {
    patch_calc_resis_nl1_fortran(p_resis);
#endif
  } else {
    assert(0);
  }
}

// ----------------------------------------------------------------------
// patch_calc_resis_const

static void _mrc_unused
patch_calc_resis_const(fld3d_t p_resis, fld3d_t p_Jcc, fld3d_t p_U, 
		       fld3d_t p_zmask)
{
  if (s_opt_mhd_calc_resis == OPT_MHD_C) {
    patch_calc_resis_const_c(p_resis, p_Jcc, p_U, p_zmask);
#if defined(HAVE_OPENGGCM_FORTRAN) && defined(MRC_FLD_AS_FLOAT_H)
  } else if (s_opt_mhd_calc_resis == OPT_MHD_FORTRAN) {
    patch_calc_resis_const_fortran(p_resis, p_Jcc, p_U, p_zmask);
#endif
  } else {
    assert(0);
  }
}

#endif
