
// ======================================================================
// nl1

// ----------------------------------------------------------------------
// patch_calc_resis_nl1_c

static void
patch_calc_resis_nl1_c(fld3d_t p_f, int m_curr)
{
  // FIXME, unnecessary
  fld3d_foreach(i,j,k, 2, 2) {
    F3S(p_f, _RESIS, i,j,k) = 0.f;
  } fld3d_foreach_end;
}

// ======================================================================
// const

// ----------------------------------------------------------------------
// patch_res1_const

static void _mrc_unused
patch_res1_const(fld3d_t p_f)
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

// ----------------------------------------------------------------------
// patch_calc_resis_const_c

static void
patch_calc_resis_const_c(fld3d_t p_f, int m_curr)
{
  patch_curbc(p_f, m_curr);
  patch_res1_const(p_f);
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
patch_calc_resis_nl1_fortran(fld3d_t p_f, int m_curr)
{
  calc_resis_nl1_F77(F(p_f, _B1X + m_curr), F(p_f, _B1Y + m_curr), F(p_f, _B1Z + m_curr),
		     F(p_f, _RESIS));
}

static void
patch_calc_resis_const_fortran(fld3d_t p_f, int m_curr)
{
  calc_resis_const_F77(F(p_f, _B1X + m_curr), F(p_f, _B1Y + m_curr), F(p_f, _B1Z + m_curr),
		       F(p_f, _CURRX), F(p_f, _CURRY), F(p_f, _CURRZ),
		       F(p_f, _TMP1), F(p_f, _TMP2), F(p_f, _TMP3),
		       F(p_f, _FLX), F(p_f, _FLY), F(p_f, _FLZ),
		       F(p_f, _ZMASK), F(p_f, _RR), F(p_f, _PP), F(p_f, _RESIS));
}

#endif

// ======================================================================
// dispatch

// ----------------------------------------------------------------------
// patch_calc_resis_nl1

static void _mrc_unused
patch_calc_resis_nl1(fld3d_t p_f, int m_curr)
{
  if (s_opt_mhd_calc_resis == OPT_MHD_C) {
    patch_calc_resis_nl1_c(p_f, m_curr);
#if defined(HAVE_OPENGGCM_FORTRAN) && defined(MRC_FLD_AS_FLOAT_H)
  } else if (s_opt_mhd_calc_resis == OPT_MHD_FORTRAN) {
    patch_calc_resis_nl1_fortran(p_f, m_curr);
#endif
  } else {
    assert(0);
  }
}

// ----------------------------------------------------------------------
// patch_calc_resis_const

static void _mrc_unused
patch_calc_resis_const(fld3d_t p_f, int m_curr)
{
  if (s_opt_mhd_calc_resis == OPT_MHD_C) {
    patch_calc_resis_const_c(p_f, m_curr);
#if defined(HAVE_OPENGGCM_FORTRAN) && defined(MRC_FLD_AS_FLOAT_H)
  } else if (s_opt_mhd_calc_resis == OPT_MHD_FORTRAN) {
    patch_calc_resis_const_fortran(p_f, m_curr);
#endif
  } else {
    assert(0);
  }
}

