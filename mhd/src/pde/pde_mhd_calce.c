
#ifndef PDE_MHD_CALCE_C
#define PDE_MHD_CALCE_C

// ======================================================================

static inline float
bcthy3f(mrc_fld_data_t s1, mrc_fld_data_t s2)
{
  const mrc_fld_data_t REPS = 1.e-10f;

  if (s1 > 0.f && mrc_fld_abs(s2) > REPS) {
    assert(!s_calce_aspect_low);
/* .if(calce_aspect_low) then */
/* .call lowmask(I, 0, 0,tl1) */
/* .call lowmask( 0,J, 0,tl2) */
/* .call lowmask( 0, 0,K,tl3) */
/* .call lowmask(I,J,K,tl4) */
/*       tt=tt*(1.0-max(tl1,tl2,tl3,tl4)) */
    return s1 / s2;
  }
  return 0.f;
}

#define ID(XX) (di*((XX)==0))
#define JD(XX) (dj*((XX)==1))
#define KD(XX) (dk*((XX)==2))

#define F3S_P(XX, p_f, m, i,j,k) F3S(p_f, m, i+ID(XX),j+JD(XX),k+KD(XX))
#define F3S_M(XX, p_f, m, i,j,k) F3S(p_f, m, i-ID(XX),j-JD(XX),k-KD(XX))

#define F3S_YYP(p_f, m, i,j,k) F3S_P(YY, p_f, m, i,j,k)
#define F3S_ZZP(p_f, m, i,j,k) F3S_P(ZZ, p_f, m, i,j,k)
#define F3S_YYM(p_f, m, i,j,k) F3S_M(YY, p_f, m, i,j,k)
#define F3S_ZZM(p_f, m, i,j,k) F3S_M(ZZ, p_f, m, i,j,k)

#define BT_YYP(p_B, d, i,j,k) BT_(p_B, d, i+ID(YY),j+JD(YY),k+KD(YY))
#define BT_ZZP(p_B, d, i,j,k) BT_(p_B, d, i+ID(ZZ),j+JD(ZZ),k+KD(ZZ))
#define BT_YYM(p_B, d, i,j,k) BT_(p_B, d, i-ID(YY),j-JD(YY),k-KD(YY))
#define BT_ZZM(p_B, d, i,j,k) BT_(p_B, d, i-ID(ZZ),j-JD(ZZ),k-KD(ZZ))

// ----------------------------------------------------------------------
// patch_calc_avg_dz_By

static inline void
patch_calc_avg_dz_By(fld3d_t p_dB, fld3d_t p_B, int XX, int YY, int ZZ)
{
  static fld3d_t p_tmp1;
  fld3d_setup_tmp_compat(&p_tmp1, 2, _TMP1);

  // d_z B_y, d_y B_z on x edges
  fld3d_foreach_stagger(i,j,k, 1, 2) {
#if OPT_STAGGER == OPT_STAGGER_GGCM
    mrc_fld_data_t bd1[3] = { PDE_INV_DXF(i+1), PDE_INV_DYF(j+1), PDE_INV_DZF(k+1) };

    F3S(p_tmp1, 0, i,j,k) = bd1[ZZ] * (BT_ZZP(p_B, YY, i,j,k) - BT_(p_B, YY, i,j,k));
    F3S(p_tmp1, 1, i,j,k) = bd1[YY] * (BT_YYP(p_B, ZZ, i,j,k) - BT_(p_B, ZZ, i,j,k));
#else
    mrc_fld_data_t bd1[3] = { PDE_INV_DXF(i), PDE_INV_DYF(j), PDE_INV_DZF(k) };

    F3S(p_tmp1, 0, i,j,k) = bd1[ZZ] * (BT_(p_B, YY, i,j,k) - BT_ZZM(p_B, YY, i,j,k));
    F3S(p_tmp1, 1, i,j,k) = bd1[YY] * (BT_(p_B, ZZ, i,j,k) - BT_YYM(p_B, ZZ, i,j,k));
#endif
  } fld3d_foreach_end;

  // .5 * harmonic average if same sign
  fld3d_foreach(i,j,k, 1, 1) {
    mrc_fld_data_t s1, s2;
#if OPT_STAGGER == OPT_STAGGER_GGCM
    // dz_By on y face
    s1 = F3S(p_tmp1, 0, i,j,k) * F3S_ZZM(p_tmp1, 0, i,j,k);
    s2 = F3S(p_tmp1, 0, i,j,k) + F3S_ZZM(p_tmp1, 0, i,j,k);
    F3S(p_dB, 0, i,j,k) = bcthy3f(s1, s2);
    // dy_Bz on z face
    s1 = F3S(p_tmp1, 1, i,j,k) * F3S_YYM(p_tmp1, 1, i,j,k);
    s2 = F3S(p_tmp1, 1, i,j,k) + F3S_YYM(p_tmp1, 1, i,j,k);
    F3S(p_dB, 1, i,j,k) = bcthy3f(s1, s2);
#else
    // dz_By on y face
    s1 = F3S_ZZP(p_tmp1, 0, i,j,k) * F3S(p_tmp1, 0, i,j,k);
    s2 = F3S_ZZP(p_tmp1, 0, i,j,k) + F3S(p_tmp1, 0, i,j,k);
    F3S(p_dB, 0, i,j,k) = bcthy3f(s1, s2);
    // dy_Bz on z face
    s1 = F3S_YYP(p_tmp1, 1, i,j,k) * F3S(p_tmp1, 1, i,j,k);
    s2 = F3S_YYP(p_tmp1, 1, i,j,k) + F3S(p_tmp1, 1, i,j,k);
    F3S(p_dB, 1, i,j,k) = bcthy3f(s1, s2);
#endif
  } fld3d_foreach_end;
}
 
// ----------------------------------------------------------------------
// calc_ve_x_B
//
// calculates (v - d_i J) x B, so without Hall it's just v x B

static inline void
calc_ve_x_B(mrc_fld_data_t ttmp[2], fld3d_t p_B, fld3d_t p_W, fld3d_t p_dB,
	   int i, int j, int k, int XX, int YY, int ZZ,
	   mrc_fld_data_t dt)
{
#if OPT_STAGGER == OPT_STAGGER_GGCM
  mrc_fld_data_t bd2[3] = { PDE_DX(i), PDE_DY(j), PDE_DZ(k) };
  mrc_fld_data_t bd2p[3] = { PDE_DX(i+1), PDE_DY(j+1), PDE_DZ(k+1) };
#else
  mrc_fld_data_t bd2m[3] = { PDE_DX(i-1), PDE_DY(j-1), PDE_DZ(k-1) };
  mrc_fld_data_t bd2[3] = { PDE_DX(i), PDE_DY(j), PDE_DZ(k) };
#endif

  // edge centered velocity
  mrc_fld_data_t vvYY = CC_TO_EC(p_W, VX + YY, i,j,k, XX);
  if (s_opt_hall != OPT_HALL_NONE) {
    mrc_fld_data_t vcurrYY = CC_TO_EC(s_p_aux.Jcc, YY, i, j, k, XX);
    vvYY -= s_d_i * vcurrYY;
  }

  mrc_fld_data_t vbZZ;
#if OPT_STAGGER == OPT_STAGGER_GGCM
  if (vvYY > 0.f) {
    vbZZ =     BT_(p_B, ZZ, i,j,k) +     F3S(p_dB, 1, i,j,k) * (bd2[YY]  - dt*vvYY);
  } else {
    vbZZ = BT_YYP(p_B, ZZ, i,j,k) - F3S_YYP(p_dB, 1, i,j,k) * (bd2p[YY] + dt*vvYY);
  }
#else
  if (vvYY > 0.f) {
    vbZZ = BT_YYM(p_B, ZZ, i,j,k) + F3S_YYM(p_dB, 1, i,j,k) * (bd2m[YY] - dt*vvYY);
  } else {
    vbZZ =     BT_(p_B, ZZ, i,j,k) -     F3S(p_dB, 1, i,j,k) * (bd2[YY]  + dt*vvYY);
  }
#endif
  ttmp[0] = vbZZ * vvYY;
  
  // edge centered velocity
  mrc_fld_data_t vvZZ = CC_TO_EC(p_W, VX + ZZ, i,j,k, XX);
  if (s_opt_hall != OPT_HALL_NONE) {
    mrc_fld_data_t vcurrZZ = CC_TO_EC(s_p_aux.Jcc, ZZ, i, j, k, XX);
    vvZZ -= s_d_i * vcurrZZ;
  }

  mrc_fld_data_t vbYY;
#if OPT_STAGGER == OPT_STAGGER_GGCM
  if (vvZZ > 0.f) {
    vbYY =     BT_(p_B, YY, i,j,k) +     F3S(p_dB, 0, i,j,k) * (bd2[ZZ]  - dt*vvZZ);
  } else {
    vbYY = BT_ZZP(p_B, YY, i,j,k) - F3S_ZZP(p_dB, 0, i,j,k) * (bd2p[ZZ] + dt*vvZZ);
  }
#else
  if (vvZZ > 0.f) {
    vbYY = BT_ZZM(p_B, YY, i,j,k) + F3S_ZZM(p_dB, 0, i,j,k) * (bd2m[ZZ] - dt*vvZZ);
  } else {
    vbYY =     BT_(p_B, YY, i,j,k) -     F3S(p_dB, 0, i,j,k) * (bd2[ZZ]  + dt*vvZZ);
  }
#endif
  ttmp[1] = vbYY * vvZZ;
}

// ----------------------------------------------------------------------
// patch_bcthy3z_NL1

static void
patch_bcthy3z_NL1(fld3d_t p_E, mrc_fld_data_t dt, fld3d_t p_U, fld3d_t p_W,
		  fld3d_t p_rmask,
		  int XX, int YY, int ZZ)
{
  const mrc_fld_data_t REPS = 1.e-10f;
  static fld3d_t p_dB;
  fld3d_setup_tmp_compat(&p_dB, 2, _TMP3);
  fld3d_t p_B = fld3d_make_view(p_U, BX);

  patch_calc_avg_dz_By(p_dB, p_B, XX, YY, ZZ);

  mrc_fld_data_t diffmul = 1.f;
  if (s_mhd_time < s_diff_timelo) { // no anomalous res at startup
    diffmul = 0.f;
  }

  // edge centered E = - v x B (+ dissipation)
  fld3d_foreach_stagger(i,j,k, 0, 1) {
    mrc_fld_data_t ttmp[2];
    calc_ve_x_B(ttmp, p_B, p_W, p_dB, i, j, k, XX, YY, ZZ, dt);

    mrc_fld_data_t t1m = BT_YYP(p_B, ZZ, i,j,k) - BT_(p_B, ZZ, i,j,k);
    mrc_fld_data_t t1p = mrc_fld_abs(BT_YYP(p_B, ZZ, i,j,k)) + mrc_fld_abs(BT_(p_B, ZZ, i,j,k));
    mrc_fld_data_t t2m = F3S_ZZP(p_B, YY, i,j,k) - F3S(p_B, YY, i,j,k);
    mrc_fld_data_t t2p = mrc_fld_abs(BT_ZZP(p_B, YY, i,j,k)) + mrc_fld_abs(BT_(p_B, YY, i,j,k));
    mrc_fld_data_t tp = t1p + t2p + REPS;
    mrc_fld_data_t tpi = diffmul / tp;
    mrc_fld_data_t d1 = sqr(t1m * tpi);
    mrc_fld_data_t d2 = sqr(t2m * tpi);
    if (d1 < s_diffth) d1 = 0.;
    if (d2 < s_diffth) d2 = 0.;
    ttmp[0] -= d1 * t1m * F3S(p_rmask, 0, i,j,k);
    ttmp[1] -= d2 * t2m * F3S(p_rmask, 0, i,j,k);
    //F3S(p_resis, 0, i,j,k) += mrc_fld_abs(d1+d2) * F3S(p_zmask, 0, i,j,k);
    F3S(p_E, XX, i,j,k) = ttmp[0] - ttmp[1];
  } fld3d_foreach_end;
}

// ----------------------------------------------------------------------
// patch_bcthy3z_const

static void
patch_bcthy3z_const(fld3d_t p_E, mrc_fld_data_t dt, fld3d_t p_U, fld3d_t p_W,
		    fld3d_t p_resis,
		    int XX, int YY, int ZZ)
{
  static fld3d_t p_dB;
  fld3d_setup_tmp_compat(&p_dB, 2, _TMP3);
  fld3d_t p_B = fld3d_make_view(p_U, BX);

  patch_calc_avg_dz_By(p_dB, p_B, XX, YY, ZZ);

  // edge centered E = - ve x B (+ dissipation)
  fld3d_foreach_stagger(i,j,k, 0, 1) {
    mrc_fld_data_t ttmp[2];
    calc_ve_x_B(ttmp, p_B, p_W, p_dB, i, j, k, XX, YY, ZZ, dt);

    mrc_fld_data_t vcurrXX = CC_TO_EC(s_p_aux.Jcc, XX, i,j,k, XX);
    mrc_fld_data_t vresis = CC_TO_EC(p_resis, 0, i,j,k, XX);
    F3S(p_E, XX, i,j,k) = ttmp[0] - ttmp[1] - vresis * vcurrXX;
  } fld3d_foreach_end;
}

// ----------------------------------------------------------------------
// patch_calce_nl1_c

static void
patch_calce_nl1_c(fld3d_t p_E, mrc_fld_data_t dt, fld3d_t p_U, fld3d_t p_W,
		  fld3d_t p_rmask)
{
  patch_bcthy3z_NL1(p_E, dt, p_U, p_W, p_rmask, 0,1,2);
  patch_bcthy3z_NL1(p_E, dt, p_U, p_W, p_rmask, 1,2,0);
  patch_bcthy3z_NL1(p_E, dt, p_U, p_W, p_rmask, 2,0,1);
}

// ----------------------------------------------------------------------
// patch_calce_const_c

static void
patch_calce_const_c(fld3d_t p_E, mrc_fld_data_t dt, fld3d_t p_U, fld3d_t p_W, fld3d_t p_resis)
{
  patch_bcthy3z_const(p_E, dt, p_U, p_W, p_resis, 0,1,2);
  patch_bcthy3z_const(p_E, dt, p_U, p_W, p_resis, 1,2,0);
  patch_bcthy3z_const(p_E, dt, p_U, p_W, p_resis, 2,0,1);
}

// ----------------------------------------------------------------------
// patch_calce_c

static void
patch_calce_c(fld3d_t p_E, mrc_fld_data_t dt, fld3d_t p_U, fld3d_t p_W,
	      fld3d_t p_zmask, fld3d_t p_rmask, fld3d_t p_resis,
	      fld3d_t p_Jcc)
{
  s_p_aux.Jcc = p_Jcc; // FIXME

  switch (s_magdiffu) {
  case MAGDIFFU_NL1:
    return patch_calce_nl1_c(p_E, dt, p_U, p_W, p_rmask);
  case MAGDIFFU_CONST:
    return patch_calce_const_c(p_E, dt, p_U, p_W, p_resis);
  default:
    assert(0);
  }
}

// ----------------------------------------------------------------------
// patch_calce_fortran

#if defined(HAVE_OPENGGCM_FORTRAN) && defined(MRC_FLD_AS_FLOAT_H)

#define calce_F77 F77_FUNC(calce,CALCE)

void calce_F77(real *bx2, real *by2, real *bz2,
	       real *zmask, real *rmask, real *resis,
	       real *flx, real *fly, real *flz,
	       real *vx, real *vy, real *vz,
	       real *currx, real *curry, real *currz,
	       real *dt, real *time);

static void
patch_calce_fortran(fld3d_t p_E, mrc_fld_data_t dt, fld3d_t p_U, fld3d_t p_W,
		    fld3d_t p_zmask, fld3d_t p_rmask, fld3d_t p_resis,
		    fld3d_t p_Jcc)
{
  calce_F77(F(p_U, BX), F(p_U, BY), F(p_U, BZ), 
	    F(p_zmask, 0), F(p_rmask, 0), F(p_resis, 0),
	    F(p_E, 0), F(p_E, 1), F(p_E, 2), 
	    F(p_W, VX), F(p_W, VY), F(p_W, VZ), 
	    F(p_Jcc, 0), F(p_Jcc, 1), F(p_Jcc, 2), 
	    &dt, &s_mhd_time);
}

#endif

// ----------------------------------------------------------------------
// patch_calce

static void _mrc_unused
patch_calce(fld3d_t p_E, mrc_fld_data_t dt, fld3d_t p_U, fld3d_t p_W,
	    fld3d_t p_zmask, fld3d_t p_rmask, fld3d_t p_resis,
	    fld3d_t p_Jcc)
{
  if (s_opt_mhd_calce == OPT_MHD_C) {
    patch_calce_c(p_E, dt, p_U, p_W, p_zmask, p_rmask, p_resis, p_Jcc);
#if defined(HAVE_OPENGGCM_FORTRAN) && defined(MRC_FLD_AS_FLOAT_H)
  } else if (s_opt_mhd_calce == OPT_MHD_FORTRAN) {
    patch_calce_fortran(p_E, dt, p_U, p_W, p_zmask, p_rmask, p_resis, p_Jcc);
#endif
  } else {
    assert(0);
  }
}

// ----------------------------------------------------------------------
// patch_calc_e
//
// alternate (newer, better) way, that incorporates all the steps necessary in
// actually calculating the E field, beginning to end

static void _mrc_unused
patch_calc_e(fld3d_t p_E, mrc_fld_data_t dt, fld3d_t p_U, fld3d_t p_W,
	     fld3d_t p_zmask, fld3d_t p_rmask)
{
  static fld3d_t p_resis;
  fld3d_setup_tmp_compat(&p_resis, 1, _RESIS);

  if (s_opt_hall != OPT_HALL_NONE || s_magdiffu == MAGDIFFU_CONST) {
    fld3d_setup_tmp(&s_p_aux.Jcc, 3);
    patch_calc_current_cc(s_p_aux.Jcc, p_U, p_zmask);
  }

  switch (s_magdiffu) {
  case MAGDIFFU_NL1:
    patch_calce_nl1_c(p_E, dt, p_U, p_W, p_rmask);
    break;
  case MAGDIFFU_RES1:
    assert(0);
    // calc_resis_res1(bxB,byB,bzB,currx,curry,currz,tmp1,tmp2,tmp3,flx,fly,flz,zmask,rr,pp,resis);
    // calce...
    break;
  case MAGDIFFU_CONST:
    patch_res1_const(p_resis);
    patch_calce_const_c(p_E, dt, p_U, p_W, p_resis);
    break;
  default:
    assert(0);
  }
}

#endif

