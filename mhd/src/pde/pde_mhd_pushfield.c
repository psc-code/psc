
// ----------------------------------------------------------------------
// patch_calc_resis_nl1_c

static void
patch_calc_resis_nl1_c(fld3d_t p_f, int m_curr)
{
  fld3d_foreach(i,j,k, 2, 2) {
    F3S(p_f, _RESIS, i,j,k) = 0.f;
  } fld3d_foreach_end;
}

// ----------------------------------------------------------------------
// patch_calc_resis_nl1_fortran

#if defined(HAVE_OPENGGCM_FORTRAN) && defined(MRC_FLD_AS_FLOAT_H)

#define calc_resis_nl1_F77 F77_FUNC(calc_resis_nl1,CALC_RESIS_NL1)

void calc_resis_nl1_F77(real *bx, real *by, real *bz, real *resis);

static void
patch_calc_resis_nl1_fortran(fld3d_t p_f, int m_curr)
{
  calc_resis_nl1_F77(F(p_f, _B1X + m_curr), F(p_f, _B1Y + m_curr), F(p_f, _B1Z + m_curr),
		     F(p_f, _RESIS));
}

#endif

// ----------------------------------------------------------------------
// patch_calc_resis_nl1

static void
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
// vgr0

static void
vgr0(fld3d_t p_f, int m)
{
  fld3d_foreach(ix,iy,iz, 2, 2) {
    F3S(p_f, m, ix,iy,iz) = 0.;
  } fld3d_foreach_end;
}

// ----------------------------------------------------------------------
// patch_pushfield1_c

static void
patch_pushfield1_c(fld3d_t p_f, mrc_fld_data_t dt)
{
  if (s_magdiffu == MAGDIFFU_NL1) {
    patch_calc_resis_nl1(p_f, _RR1);
    vgr0(p_f, _CURRX);
    vgr0(p_f, _CURRY);
    vgr0(p_f, _CURRZ);
  } else if (s_magdiffu == MAGDIFFU_RES1) {
    assert(0);
    //    calc_resis_res1(bxB,byB,bzB,currx,curry,currz,tmp1,tmp2,tmp3,flx,fly,flz,zmask,rr,pp,resis);
  } else if (s_magdiffu == MAGDIFFU_CONST) {
    assert(0);
    //calc_resis_const(bxB,byB,bzB,currx,curry,currz,tmp1,tmp2,tmp3,flx,fly,flz,zmask,rr,pp,resis);
  }

  patch_push_ej(p_f, dt, _RR1, _RR2);
  patch_pfie3(p_f, dt, _RR1, _RR1, _RR2);
}

// ----------------------------------------------------------------------
// patch_pushfield1_fortran

#if defined(HAVE_OPENGGCM_FORTRAN) && defined(MRC_FLD_AS_FLOAT_H)

#define pushfield1_F77 F77_FUNC(pushfield1,PUSHFIELD1)

void pushfield1_F77(real *rr1, real *rv1x, real *rv1y, real *rv1z, real *uu1,
		    real *b1x, real *b1y, real *b1z,
		    real *rr2, real *rv2x, real *rv2y, real *rv2z, real *uu2,
		    real *b2x, real *b2y, real *b2z,
		    real *rr, real *vx, real *vy, real *vz, real *pp,
		    real *cmsv, real *ymask, real *zmask, real *rmask,
		    real *flx, real *fly, real *flz,
		    real *tmp1, real *tmp2, real *tmp3, real *resis,
		    real *dth, real *time);

static void
patch_pushfield1_fortran(fld3d_t p_f, mrc_fld_data_t dt)
{
  pushfield1_F77(F(p_f, _RR1), F(p_f, _RV1X), F(p_f, _RV1Y), F(p_f, _RV1Z), F(p_f, _UU1),
		 F(p_f, _B1X), F(p_f, _B1Y), F(p_f, _B1Z),
		 F(p_f, _RR2), F(p_f, _RV2X), F(p_f, _RV2Y), F(p_f, _RV2Z), F(p_f, _UU2),
		 F(p_f, _B2X), F(p_f, _B2Y), F(p_f, _B2Z),
		 F(p_f, _RR), F(p_f, _VX), F(p_f, _VY), F(p_f, _VZ), F(p_f, _PP),
		 F(p_f, _CMSV), F(p_f, _YMASK), F(p_f, _ZMASK), F(p_f, _RMASK),
		 F(p_f, _FLX), F(p_f, _FLY), F(p_f, _FLZ),
		 F(p_f, _TMP1), F(p_f, _TMP2), F(p_f, _TMP3), F(p_f, _RESIS),
		 &dt, &s_mhd_time);
}

#endif

// ----------------------------------------------------------------------
// patch_pushfield1

static void _mrc_unused
patch_pushfield1(fld3d_t p_f, mrc_fld_data_t dt)
{
  if (s_opt_mhd_pushfield1 == OPT_MHD_C) {
    patch_pushfield1_c(p_f, dt);
#if defined(HAVE_OPENGGCM_FORTRAN) && defined(MRC_FLD_AS_FLOAT_H)
  } else if (s_opt_mhd_pushfield1 == OPT_MHD_FORTRAN) {
    patch_pushfield1_fortran(p_f, dt);
#endif
  } else {
    assert(0);
  }
}

// ======================================================================

// ----------------------------------------------------------------------
// patch_pushfield2_c

static void
patch_pushfield2_c(fld3d_t p_f, mrc_fld_data_t dt)
{
  if (s_magdiffu == MAGDIFFU_NL1) {
    patch_calc_resis_nl1(p_f, _RR2);
    vgr0(p_f, _CURRX);
    vgr0(p_f, _CURRY);
    vgr0(p_f, _CURRZ);
  } else if (s_magdiffu == MAGDIFFU_RES1) {
    assert(0);
    //    calc_resis_res1(bxB,byB,bzB,currx,curry,currz,tmp1,tmp2,tmp3,flx,fly,flz,zmask,rr,pp,resis);
  } else if (s_magdiffu == MAGDIFFU_CONST) {
    assert(0);
    //calc_resis_const(bxB,byB,bzB,currx,curry,currz,tmp1,tmp2,tmp3,flx,fly,flz,zmask,rr,pp,resis);
  }

  patch_push_ej(p_f, dt, _RR2, _RR1);
  patch_pfie3(p_f, dt, _RR1, _RR2, _RR1);
}

// ----------------------------------------------------------------------
// patch_pushfield2_fortran

#if defined(HAVE_OPENGGCM_FORTRAN) && defined(MRC_FLD_AS_FLOAT_H)

#define pushfield2_F77 F77_FUNC(pushfield2,PUSHFIELD2)

void pushfield2_F77(real *rr1, real *rv1x, real *rv1y, real *rv1z, real *uu1,
		    real *b1x, real *b1y, real *b1z,
		    real *rr2, real *rv2x, real *rv2y, real *rv2z, real *uu2,
		    real *b2x, real *b2y, real *b2z,
		    real *rr, real *vx, real *vy, real *vz, real *pp,
		    real *cmsv, real *ymask, real *zmask, real *rmask,
		    real *flx, real *fly, real *flz,
		    real *tmp1, real *tmp2, real *tmp3, real *resis,
		    real *dth, real *time);

static void
patch_pushfield2_fortran(fld3d_t p_f, mrc_fld_data_t dt)
{
  pushfield2_F77(F(p_f, _RR1), F(p_f, _RV1X), F(p_f, _RV1Y), F(p_f, _RV1Z), F(p_f, _UU1),
		 F(p_f, _B1X), F(p_f, _B1Y), F(p_f, _B1Z),
		 F(p_f, _RR2), F(p_f, _RV2X), F(p_f, _RV2Y), F(p_f, _RV2Z), F(p_f, _UU2),
		 F(p_f, _B2X), F(p_f, _B2Y), F(p_f, _B2Z),
		 F(p_f, _RR), F(p_f, _VX), F(p_f, _VY), F(p_f, _VZ), F(p_f, _PP),
		 F(p_f, _CMSV), F(p_f, _YMASK), F(p_f, _ZMASK), F(p_f, _RMASK),
		 F(p_f, _FLX), F(p_f, _FLY), F(p_f, _FLZ),
		 F(p_f, _TMP1), F(p_f, _TMP2), F(p_f, _TMP3), F(p_f, _RESIS),
		 &dt, &s_mhd_time);
}

#endif

// ----------------------------------------------------------------------
// patch_pushfield2

static void _mrc_unused
patch_pushfield2(fld3d_t p_f, mrc_fld_data_t dt)
{
  if (s_opt_mhd_pushfield2 == OPT_MHD_C) {
    patch_pushfield2_c(p_f, dt);
#if defined(HAVE_OPENGGCM_FORTRAN) && defined(MRC_FLD_AS_FLOAT_H)
  } else if (s_opt_mhd_pushfield2 == OPT_MHD_FORTRAN) {
    patch_pushfield2_fortran(p_f, dt);
#endif
  } else {
    assert(0);
  }
}

