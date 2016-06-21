
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

#if OPT_STAGGER == OPT_STAGGER_GGCM

static inline void
calc_avg_dz_By(fld3d_t p_dB, fld3d_t p_U, int XX, int YY, int ZZ,
	       int JX1, int JY1, int JZ1, int JX2, int JY2, int JZ2)
{
  fld3d_t p_tmp1 = fld3d_make_tmp(2, _TMP1);

  // d_z B_y, d_y B_z on x edges
  fld3d_foreach(i,j,k, 2, 1) {
    mrc_fld_data_t bd1[3] = { PDE_INV_DXF(i+1), PDE_INV_DYF(j+1), PDE_INV_DZF(k+1) };

    F3S(p_tmp1, 0, i,j,k) = bd1[ZZ] * 
      (F3S(p_U, BX + YY, i+JX2,j+JY2,k+JZ2) - F3S(p_U, BX + YY, i,j,k));
    F3S(p_tmp1, 1, i,j,k) = bd1[YY] * 
      (F3S(p_U, BX + ZZ, i+JX1,j+JY1,k+JZ1) - F3S(p_U, BX + ZZ, i,j,k));
  } fld3d_foreach_end;

  // .5 * harmonic average if same sign
  fld3d_foreach(i,j,k, 1, 1) {
    mrc_fld_data_t s1, s2;
    // dz_By on y face
    s1 = F3S(p_tmp1, 0, i,j,k) * F3S(p_tmp1, 0, i-JX2,j-JY2,k-JZ2);
    s2 = F3S(p_tmp1, 0, i,j,k) + F3S(p_tmp1, 0, i-JX2,j-JY2,k-JZ2);
    F3S(p_dB, 0, i,j,k) = bcthy3f(s1, s2);
    // dy_Bz on z face
    s1 = F3S(p_tmp1, 1, i,j,k) * F3S(p_tmp1, 1, i-JX1,j-JY1,k-JZ1);
    s2 = F3S(p_tmp1, 1, i,j,k) + F3S(p_tmp1, 1, i-JX1,j-JY1,k-JZ1);
    F3S(p_dB, 1, i,j,k) = bcthy3f(s1, s2);
  } fld3d_foreach_end;
}

#else

static inline void
calc_avg_dz_By(fld3d_t p_dB, fld3d_t p_U, int XX, int YY, int ZZ,
	       int JX1, int JY1, int JZ1, int JX2, int JY2, int JZ2)
{
  fld3d_t p_tmp1 = fld3d_make_tmp(2, _TMP1);

  // d_z B_y, d_y B_z on x edges
  fld3d_foreach(i,j,k, 1, 2) {
    mrc_fld_data_t bd1[3] = { PDE_INV_DXF(i), PDE_INV_DYF(j), PDE_INV_DZF(k) };

    F3S(p_tmp1, 0, i,j,k) = bd1[ZZ] * 
      (F3S(p_U, BX + YY, i,j,k) - F3S(p_U, BX + YY, i-JX2,j-JY2,k-JZ2));
    F3S(p_tmp1, 1, i,j,k) = bd1[YY] * 
      (F3S(p_U, BX + ZZ, i,j,k) - F3S(p_U, BX + ZZ, i-JX1,j-JY1,k-JZ1));
  } fld3d_foreach_end;

  // .5 * harmonic average if same sign
  fld3d_foreach(i,j,k, 1, 1) {
    mrc_fld_data_t s1, s2;
    // dz_By on y face
    s1 = F3S(p_tmp1, 0, i+JX2,j+JY2,k+JZ2) * F3S(p_tmp1, 0, i,j,k);
    s2 = F3S(p_tmp1, 0, i+JX2,j+JY2,k+JZ2) + F3S(p_tmp1, 0, i,j,k);
    F3S(p_dB, 0, i,j,k) = bcthy3f(s1, s2);
    // dy_Bz on z face
    s1 = F3S(p_tmp1, 1, i+JX1,j+JY1,k+JZ1) * F3S(p_tmp1, 1, i,j,k);
    s2 = F3S(p_tmp1, 1, i+JX1,j+JY1,k+JZ1) + F3S(p_tmp1, 1, i,j,k);
    F3S(p_dB, 1, i,j,k) = bcthy3f(s1, s2);
  } fld3d_foreach_end;
}

#endif

#if OPT_STAGGER == OPT_STAGGER_GGCM

#define CC_TO_EC(p_f, m, i,j,k, I,J,K) \
  (.25f * (F3S(p_f, m, i  ,j  ,k  ) +  \
	   F3S(p_f, m, i  ,j+J,k+K) +  \
	   F3S(p_f, m, i+I,j  ,k+K) +  \
	   F3S(p_f, m, i+I,j+J,k  )))

#else

#define CC_TO_EC(f, m, i,j,k, I,J,K) \
  (.25f * (F3(f, m, i-I,j-J,k-K) +  \
	   F3(f, m, i-I,j  ,k  ) +  \
	   F3(f, m, i  ,j-J,k  ) +  \
	   F3(f, m, i  ,j  ,k-K)))

#endif

#if OPT_STAGGER == OPT_STAGGER_GGCM

static inline void
calc_v_x_B(mrc_fld_data_t ttmp[2], fld3d_t p_U, fld3d_t p_W, fld3d_t p_dB,
	   int i, int j, int k,
	   int XX, int YY, int ZZ, int I, int J, int K,
	   int JX1, int JY1, int JZ1, int JX2, int JY2, int JZ2,
	   mrc_fld_data_t dt)
{
  mrc_fld_data_t bd2[3] = { PDE_DX(i), PDE_DY(j), PDE_DZ(k) };
  mrc_fld_data_t bd2p[3] = { PDE_DX(i+1), PDE_DY(j+1), PDE_DZ(k+1) };
  mrc_fld_data_t vbZZ;
  // edge centered velocity
  mrc_fld_data_t vvYY = CC_TO_EC(p_W, VX + YY, i,j,k, I,J,K) /* - d_i * vcurrYY */;
  if (vvYY > 0.f) {
    vbZZ = F3S(p_U, BX + ZZ, i,j,k) +
      F3S(p_dB, 1, i,j,k) * (bd2[YY] - dt*vvYY);
  } else {
    vbZZ = F3S(p_U, BX + ZZ, i+JX1,j+JY1,k+JZ1) -
      F3S(p_dB, 1, i+JX1,j+JY1,k+JZ1) * (bd2p[YY] + dt*vvYY);
  }
  ttmp[0] = vbZZ * vvYY;
  
  mrc_fld_data_t vbYY;
  // edge centered velocity
  mrc_fld_data_t vvZZ = CC_TO_EC(p_W, VX + ZZ, i,j,k, I,J,K) /* - d_i * vcurrZZ */;
  if (vvZZ > 0.f) {
    vbYY = F3S(p_U, BX + YY, i,j,k) +
      F3S(p_dB, 0, i,j,k) * (bd2[ZZ] - dt*vvZZ);
  } else {
    vbYY = F3S(p_U, BX + YY, i+JX2,j+JY2,k+JZ2) -
      F3S(p_dB, 0, i+JX2,j+JY2,k+JZ2) * (bd2p[ZZ] + dt*vvZZ);
  }
  ttmp[1] = vbYY * vvZZ;
}

#else

static inline void
calc_v_x_B(mrc_fld_data_t ttmp[2], struct mrc_fld *f, int m_curr, int i, int j, int k,
	   int XX, int YY, int ZZ, int I, int J, int K,
	   int JX1, int JY1, int JZ1, int JX2, int JY2, int JZ2,
	   float *bd2x, float *bd2y, float *bd2z, mrc_fld_data_t dt)
{
    mrc_fld_data_t bd2m[3] = { bd2x[i-1], bd2y[j-1], bd2z[k-1] };
    mrc_fld_data_t bd2[3] = { bd2x[i], bd2y[j], bd2z[k] };
    mrc_fld_data_t vbZZ;
    // edge centered velocity
    mrc_fld_data_t vvYY = CC_TO_EC(f, _VX + YY, i,j,k, I,J,K) /* - d_i * vcurrYY */;
    if (vvYY > 0.f) {
      vbZZ = F3(f, m_curr + _B1X + ZZ, i-JX1,j-JY1,k-JZ1) +
	F3(f, _TMP4, i-JX1,j-JY1,k-JZ1) * (bd2m[YY] - dt*vvYY);
    } else {
      vbZZ = F3(f, m_curr + _B1X + ZZ, i,j,k) -
	F3(f, _TMP4, i,j,k) * (bd2[YY] + dt*vvYY);
    }
    ttmp[0] = vbZZ * vvYY;

    mrc_fld_data_t vbYY;
    // edge centered velocity
    mrc_fld_data_t vvZZ = CC_TO_EC(f, _VX + ZZ, i,j,k, I,J,K) /* - d_i * vcurrZZ */;
    if (vvZZ > 0.f) {
      vbYY = F3(f, m_curr + _B1X + YY, i-JX2,j-JY2,k-JZ2) +
	F3(f, _TMP3, i-JX2,j-JY2,k-JZ2) * (bd2m[ZZ] - dt*vvZZ);
    } else {
      vbYY = F3(f, m_curr + _B1X + YY, i,j,k) -
	F3(f, _TMP3, i,j,k) * (bd2[ZZ] + dt*vvZZ);
    }
    ttmp[1] = vbYY * vvZZ;

}

#endif

#if OPT_STAGGER == OPT_STAGGER_GGCM

static void
bcthy3z_NL1(fld3d_t p_E, mrc_fld_data_t dt, fld3d_t p_U, fld3d_t p_W, fld3d_t p_zmask,
	    fld3d_t p_rmask, fld3d_t p_resis, int XX, int YY, int ZZ, int I, int J, int K,
	    int JX1, int JY1, int JZ1, int JX2, int JY2, int JZ2)
{
  const mrc_fld_data_t REPS = 1.e-10f;
  fld3d_t p_dB = fld3d_make_tmp(2, _TMP3);

  calc_avg_dz_By(p_dB, p_U, XX, YY, ZZ, JX1, JY1, JZ1, JX2, JY2, JZ2);

  mrc_fld_data_t diffmul = 1.f;
  if (s_mhd_time < s_diff_timelo) { // no anomalous res at startup
    diffmul = 0.f;
  }

  // edge centered E = - v x B (+ dissipation)
  fld3d_foreach(i,j,k, 1, 0) {
    mrc_fld_data_t ttmp[2];
    calc_v_x_B(ttmp, p_U, p_W, p_dB, i, j, k, XX, YY, ZZ, I, J, K,
	       JX1, JY1, JZ1, JX2, JY2, JZ2, dt);
    
    mrc_fld_data_t t1m = F3S(p_U, BX + ZZ, i+JX1,j+JY1,k+JZ1) - F3S(p_U, BX + ZZ, i,j,k);
    mrc_fld_data_t t1p = mrc_fld_abs(F3S(p_U, BX + ZZ, i+JX1,j+JY1,k+JZ1)) + mrc_fld_abs(F3S(p_U, BX + ZZ, i,j,k));
    mrc_fld_data_t t2m = F3S(p_U, BX + YY, i+JX2,j+JY2,k+JZ2) - F3S(p_U, BX + YY, i,j,k);
    mrc_fld_data_t t2p = mrc_fld_abs(F3S(p_U, BX + YY, i+JX2,j+JY2,k+JZ2)) + mrc_fld_abs(F3S(p_U, BX + YY, i,j,k));
    mrc_fld_data_t tp = t1p + t2p + REPS;
    mrc_fld_data_t tpi = diffmul / tp;
    mrc_fld_data_t d1 = sqr(t1m * tpi);
    mrc_fld_data_t d2 = sqr(t2m * tpi);
    if (d1 < s_diffth) d1 = 0.;
    if (d2 < s_diffth) d2 = 0.;
    ttmp[0] -= d1 * t1m * F3S(p_rmask, 0, i,j,k);
    ttmp[1] -= d2 * t2m * F3S(p_rmask, 0, i,j,k);
    F3S(p_resis, 0, i,j,k) += fabsf(d1+d2) * F3S(p_zmask, 0, i,j,k);
    F3S(p_E, XX, i,j,k) = ttmp[0] - ttmp[1];
  } fld3d_foreach_end;
}

#else
static void
bcthy3z_NL1(struct ggcm_mhd *mhd, int XX, int YY, int ZZ, int I, int J, int K,
	    int JX1, int JY1, int JZ1, int JX2, int JY2, int JZ2,
	    mrc_fld_data_t dt, int m_curr)
{
  struct mrc_fld *f = mhd->fld;

  float *bd2x = ggcm_mhd_crds_get_crd(mhd->crds, 0, BD2);
  float *bd2y = ggcm_mhd_crds_get_crd(mhd->crds, 1, BD2);
  float *bd2z = ggcm_mhd_crds_get_crd(mhd->crds, 2, BD2);

  fld3d_t p_dB = fld3d_make_view(s_p_f, _TMP3), p_U = fld3d_make_view(s_p_f, m_curr);
  calc_avg_dz_By(p_dB, p_U, XX, YY, ZZ, JX1, JY1, JZ1, JX2, JY2, JZ2);

  mrc_fld_data_t diffmul=1.0;
  if (mhd->time < mhd->par.diff_timelo) { // no anomalous res at startup
    diffmul = 0.f;
  }

#define REPS (1.e-10f)

  // edge centered E = - v x B (+ dissipation)
  fld3d_foreach(i,j,k, 0, 1) {
    mrc_fld_data_t ttmp[2];
    calc_v_x_B(ttmp, f, m_curr, i, j, k, XX, YY, ZZ, I, J, K,
	       JX1, JY1, JZ1, JX2, JY2, JZ2, bd2x, bd2y, bd2z, dt);
    
    mrc_fld_data_t t1m = F3(f, m_curr + _B1X + ZZ, i+JX1,j+JY1,k+JZ1) - F3(f, m_curr + _B1X + ZZ, i,j,k);
    mrc_fld_data_t t1p = fabsf(F3(f, m_curr + _B1X + ZZ, i+JX1,j+JY1,k+JZ1)) + fabsf(F3(f, m_curr + _B1X + ZZ, i,j,k));
    mrc_fld_data_t t2m = F3(f, m_curr + _B1X + YY, i+JX2,j+JY2,k+JZ2) - F3(f, m_curr + _B1X + YY, i,j,k);
    mrc_fld_data_t t2p = fabsf(F3(f, m_curr + _B1X + YY, i+JX2,j+JY2,k+JZ2)) + fabsf(F3(f, m_curr + _B1X + YY, i,j,k));
    mrc_fld_data_t tp = t1p + t2p + REPS;
    mrc_fld_data_t tpi = diffmul / tp;
    mrc_fld_data_t d1 = sqr(t1m * tpi);
    mrc_fld_data_t d2 = sqr(t2m * tpi);
    if (d1 < mhd->par.diffth) d1 = 0.;
    if (d2 < mhd->par.diffth) d2 = 0.;
    ttmp[0] -= d1 * t1m * F3(f, _RMASK, i,j,k);
    ttmp[1] -= d2 * t2m * F3(f, _RMASK, i,j,k);
    F3(f, _RESIS, i,j,k) += fabsf(d1+d2) * F3(f, _ZMASK, i,j,k);
    F3(f, _FLX + XX, i,j,k) = ttmp[0] - ttmp[1];
  } fld3d_foreach_end;
}

#endif

#if OPT_STAGGER == OPT_STAGGER_GGCM

static void
bcthy3z_const(fld3d_t p_E, mrc_fld_data_t dt, fld3d_t p_U, fld3d_t p_W, fld3d_t p_resis, fld3d_t p_Jcc,
	      int XX, int YY, int ZZ, int I, int J, int K,
	      int JX1, int JY1, int JZ1, int JX2, int JY2, int JZ2)
{
  fld3d_t p_dB = fld3d_make_tmp(2, _TMP3);

  calc_avg_dz_By(p_dB, p_U, XX, YY, ZZ, JX1, JY1, JZ1, JX2, JY2, JZ2);

  // edge centered E = - v x B (+ dissipation)
  fld3d_foreach(i,j,k, 1, 0) {
    mrc_fld_data_t ttmp[2];
    calc_v_x_B(ttmp, p_U, p_W, p_dB, i, j, k, XX, YY, ZZ, I, J, K,
	       JX1, JY1, JZ1, JX2, JY2, JZ2, dt);

    mrc_fld_data_t vcurrXX = CC_TO_EC(p_Jcc, XX, i,j,k, I,J,K);
    mrc_fld_data_t vresis = CC_TO_EC(p_resis, 0, i,j,k, I,J,K);
    F3S(p_E, XX, i,j,k) = ttmp[0] - ttmp[1] - vresis * vcurrXX;
  } fld3d_foreach_end;
}

#else

static void
bcthy3z_const(struct ggcm_mhd *mhd, int XX, int YY, int ZZ, int I, int J, int K,
	      int JX1, int JY1, int JZ1, int JX2, int JY2, int JZ2, mrc_fld_data_t dt, int m_curr)
{
  struct mrc_fld *f = mhd->fld;

  float *bd2x = ggcm_mhd_crds_get_crd(mhd->crds, 0, BD2);
  float *bd2y = ggcm_mhd_crds_get_crd(mhd->crds, 1, BD2);
  float *bd2z = ggcm_mhd_crds_get_crd(mhd->crds, 2, BD2);

  fld3d_t p_dB = fld3d_make_view(s_p_f, _TMP3), p_U = fld3d_make_view(s_p_f, m_curr);
  calc_avg_dz_By(p_dB, p_U, XX, YY, ZZ, JX1, JY1, JZ1, JX2, JY2, JZ2);

  // edge centered E = - v x B (+ dissipation)
  fld3d_foreach(i,j,k, 0, 1) {
    mrc_fld_data_t ttmp[2];
    calc_v_x_B(ttmp, f, m_curr, i, j, k, XX, YY, ZZ, I, J, K,
	       JX1, JY1, JZ1, JX2, JY2, JZ2, bd2x, bd2y, bd2z, dt);

    mrc_fld_data_t vcurrXX = CC_TO_EC(f, _CURRX + XX, i,j,k, I,J,K);
    mrc_fld_data_t vresis = CC_TO_EC(f, _RESIS, i,j,k, I,J,K);
    F3(f, _FLX + XX, i,j,k) = ttmp[0] - ttmp[1] - vresis * vcurrXX;
  } fld3d_foreach_end;
}

#endif

#if OPT_STAGGER == OPT_STAGGER_GGCM

static void
calce_nl1_c(fld3d_t p_E, mrc_fld_data_t dt, fld3d_t p_U, fld3d_t p_W, fld3d_t p_zmask,
	    fld3d_t p_rmask, fld3d_t p_resis)
{
  bcthy3z_NL1(p_E, dt, p_U, p_W, p_zmask, p_rmask, p_resis, 0,1,2, 0,1,1, 0,1,0, 0,0,1);
  bcthy3z_NL1(p_E, dt, p_U, p_W, p_zmask, p_rmask, p_resis, 1,2,0, 1,0,1, 0,0,1, 1,0,0);
  bcthy3z_NL1(p_E, dt, p_U, p_W, p_zmask, p_rmask, p_resis, 2,0,1, 1,1,0, 1,0,0, 0,1,0);
}

#else
static void
calce_nl1_c(struct ggcm_mhd *mhd, mrc_fld_data_t dt, int m_curr)
{
  bcthy3z_NL1(mhd, 0,1,2, 0,1,1, 0,1,0, 0,0,1, dt, m_curr);
  bcthy3z_NL1(mhd, 1,2,0, 1,0,1, 0,0,1, 1,0,0, dt, m_curr);
  bcthy3z_NL1(mhd, 2,0,1, 1,1,0, 1,0,0, 0,1,0, dt, m_curr);
}

#endif

#if OPT_STAGGER == OPT_STAGGER_GGCM

static void
calce_const_c(fld3d_t p_E, mrc_fld_data_t dt, fld3d_t p_U, fld3d_t p_W, fld3d_t p_resis,
	      fld3d_t p_Jcc)
{
  bcthy3z_const(p_E, dt, p_U, p_W, p_resis, p_Jcc, 0,1,2, 0,1,1, 0,1,0, 0,0,1);
  bcthy3z_const(p_E, dt, p_U, p_W, p_resis, p_Jcc, 1,2,0, 1,0,1, 0,0,1, 1,0,0);
  bcthy3z_const(p_E, dt, p_U, p_W, p_resis, p_Jcc, 2,0,1, 1,1,0, 1,0,0, 0,1,0);
}

#else
static void
calce_const_c(struct ggcm_mhd *mhd, mrc_fld_data_t dt, int m_curr)
{
  bcthy3z_const(mhd, 0,1,2, 0,1,1, 0,1,0, 0,0,1, dt, m_curr);
  bcthy3z_const(mhd, 1,2,0, 1,0,1, 0,0,1, 1,0,0, dt, m_curr);
  bcthy3z_const(mhd, 2,0,1, 1,1,0, 1,0,0, 0,1,0, dt, m_curr);
}

#endif

// ----------------------------------------------------------------------
// patch_calce_c

static void
patch_calce_c(fld3d_t p_E, mrc_fld_data_t dt, fld3d_t p_U, fld3d_t p_W,
	      fld3d_t p_zmask, fld3d_t p_rmask, fld3d_t p_resis,
	      fld3d_t p_Jcc)
{
  switch (s_magdiffu) {
#if OPT_STAGGER == OPT_STAGGER_GGCM
  case MAGDIFFU_NL1:
    return calce_nl1_c(p_E, dt, p_U, p_W, p_zmask, p_rmask, p_resis);
  case MAGDIFFU_CONST:
    return calce_const_c(p_E, dt, p_U, p_W, p_resis, p_Jcc);
#endif
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

static void
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

#endif
