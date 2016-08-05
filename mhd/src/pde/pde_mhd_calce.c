
#ifndef PDE_MHD_CALCE_C
#define PDE_MHD_CALCE_C

// ======================================================================

static inline float
bcthy3f(mrc_fld_data_t s1, mrc_fld_data_t s2)
{
  const mrc_fld_data_t REPS = 1.e-10f;

  if (s1 > 0.f && fabsf(s2) > REPS) {
    assert(!s_calce_aspect_low);
/* .if(calce_aspect_low) then */
/* .call lowmask(IX, 0, 0,tl1) */
/* .call lowmask( 0,IY, 0,tl2) */
/* .call lowmask( 0, 0,IZ,tl3) */
/* .call lowmask(IX,IY,IZ,tl4) */
/*       tt=tt*(1.0-max(tl1,tl2,tl3,tl4)) */
    return s1 / s2;
  }
  return 0.f;
}

static inline void
calc_avg_dz_By(fld3d_t p_f, int m_curr, int XX, int YY, int ZZ,
	       int JX1, int JY1, int JZ1, int JX2, int JY2, int JZ2)
{
  // d_z B_y, d_y B_z on x edges
  fld3d_foreach(ix,iy,iz, 2, 1) {
    mrc_fld_data_t bd1[3] = { BD1X(ix), BD1Y(iy), BD1Z(iz) };

    F3S(p_f, _TMP1, ix,iy,iz) = bd1[ZZ] * 
      (F3S(p_f, m_curr + _B1X + YY, ix+JX2,iy+JY2,iz+JZ2) - F3S(p_f, m_curr + _B1X + YY, ix,iy,iz));
    F3S(p_f, _TMP2, ix,iy,iz) = bd1[YY] * 
      (F3S(p_f, m_curr + _B1X + ZZ, ix+JX1,iy+JY1,iz+JZ1) - F3S(p_f, m_curr + _B1X + ZZ, ix,iy,iz));
  } fld3d_foreach_end;

  // .5 * harmonic average if same sign
  fld3d_foreach(ix,iy,iz, 1, 1) {
    mrc_fld_data_t s1, s2;
    // dz_By on y face
    s1 = F3S(p_f, _TMP1, ix,iy,iz) * F3S(p_f, _TMP1, ix-JX2,iy-JY2,iz-JZ2);
    s2 = F3S(p_f, _TMP1, ix,iy,iz) + F3S(p_f, _TMP1, ix-JX2,iy-JY2,iz-JZ2);
    F3S(p_f, _TMP3, ix,iy,iz) = bcthy3f(s1, s2);
    // dy_Bz on z face
    s1 = F3S(p_f, _TMP2, ix,iy,iz) * F3S(p_f, _TMP2, ix-JX1,iy-JY1,iz-JZ1);
    s2 = F3S(p_f, _TMP2, ix,iy,iz) + F3S(p_f, _TMP2, ix-JX1,iy-JY1,iz-JZ1);
    F3S(p_f, _TMP4, ix,iy,iz) = bcthy3f(s1, s2);
  } fld3d_foreach_end;
}

#define CC_TO_EC(p_f, m, ix,iy,iz, IX,IY,IZ) \
  (.25f * (F3S(p_f, m, ix   ,iy   ,iz   ) +  \
	   F3S(p_f, m, ix   ,iy+IY,iz+IZ) +  \
	   F3S(p_f, m, ix+IX,iy   ,iz+IZ) +  \
	   F3S(p_f, m, ix+IX,iy+IY,iz   )))

static inline void
calc_v_x_B(mrc_fld_data_t ttmp[2], fld3d_t p_f, int m_curr, int ix, int iy, int iz,
	   int XX, int YY, int ZZ, int IX, int IY, int IZ,
	   int JX1, int JY1, int JZ1, int JX2, int JY2, int JZ2,
	   mrc_fld_data_t dt)
{
  mrc_fld_data_t bd2[3] = { BD2X(ix), BD2Y(iy), BD2Z(iz) };
  mrc_fld_data_t bd2p[3] = { BD2X(ix+1), BD2Y(iy+1), BD2Z(iz+1) };
  mrc_fld_data_t vbZZ;
  // edge centered velocity
  mrc_fld_data_t vvYY = CC_TO_EC(p_f, _VX + YY, ix,iy,iz, IX,IY,IZ) /* - d_i * vcurrYY */;
  if (vvYY > 0.f) {
    vbZZ = F3S(p_f, m_curr + _B1X + ZZ, ix,iy,iz) +
      F3S(p_f, _TMP4, ix,iy,iz) * (bd2[YY] - dt*vvYY);
  } else {
    vbZZ = F3S(p_f, m_curr + _B1X + ZZ, ix+JX1,iy+JY1,iz+JZ1) -
      F3S(p_f, _TMP4, ix+JX1,iy+JY1,iz+JZ1) * (bd2p[YY] + dt*vvYY);
  }
  ttmp[0] = vbZZ * vvYY;
  
  mrc_fld_data_t vbYY;
  // edge centered velocity
  mrc_fld_data_t vvZZ = CC_TO_EC(p_f, _VX + ZZ, ix,iy,iz, IX,IY,IZ) /* - d_i * vcurrZZ */;
  if (vvZZ > 0.f) {
    vbYY = F3S(p_f, m_curr + _B1X + YY, ix,iy,iz) +
      F3S(p_f, _TMP3, ix,iy,iz) * (bd2[ZZ] - dt*vvZZ);
  } else {
    vbYY = F3S(p_f, m_curr + _B1X + YY, ix+JX2,iy+JY2,iz+JZ2) -
      F3S(p_f, _TMP3, ix+JX2,iy+JY2,iz+JZ2) * (bd2p[ZZ] + dt*vvZZ);
  }
  ttmp[1] = vbYY * vvZZ;
}

static void
bcthy3z_NL1(fld3d_t p_f, int XX, int YY, int ZZ, int IX, int IY, int IZ,
	    int JX1, int JY1, int JZ1, int JX2, int JY2, int JZ2,
	    mrc_fld_data_t dt, int m_curr)
{
  const mrc_fld_data_t REPS = 1.e-10f;

  calc_avg_dz_By(p_f, m_curr, XX, YY, ZZ, JX1, JY1, JZ1, JX2, JY2, JZ2);

  mrc_fld_data_t diffmul = 1.f;
  if (s_mhd_time < s_diff_timelo) { // no anomalous res at startup
    diffmul = 0.f;
  }

  // edge centered E = - v x B (+ dissipation)
  fld3d_foreach(ix,iy,iz, 1, 0) {
    mrc_fld_data_t ttmp[2];
    calc_v_x_B(ttmp, p_f, m_curr, ix, iy, iz, XX, YY, ZZ, IX, IY, IZ,
	       JX1, JY1, JZ1, JX2, JY2, JZ2, dt);
    
    mrc_fld_data_t t1m = F3S(p_f, m_curr + _B1X + ZZ, ix+JX1,iy+JY1,iz+JZ1) - F3S(p_f, m_curr + _B1X + ZZ, ix,iy,iz);
    mrc_fld_data_t t1p = fabsf(F3S(p_f, m_curr + _B1X + ZZ, ix+JX1,iy+JY1,iz+JZ1)) + fabsf(F3S(p_f, m_curr + _B1X + ZZ, ix,iy,iz));
    mrc_fld_data_t t2m = F3S(p_f, m_curr + _B1X + YY, ix+JX2,iy+JY2,iz+JZ2) - F3S(p_f, m_curr + _B1X + YY, ix,iy,iz);
    mrc_fld_data_t t2p = fabsf(F3S(p_f, m_curr + _B1X + YY, ix+JX2,iy+JY2,iz+JZ2)) + fabsf(F3S(p_f, m_curr + _B1X + YY, ix,iy,iz));
    mrc_fld_data_t tp = t1p + t2p + REPS;
    mrc_fld_data_t tpi = diffmul / tp;
    mrc_fld_data_t d1 = sqr(t1m * tpi);
    mrc_fld_data_t d2 = sqr(t2m * tpi);
    if (d1 < s_diffth) d1 = 0.;
    if (d2 < s_diffth) d2 = 0.;
    ttmp[0] -= d1 * t1m * F3S(p_f, _RMASK, ix,iy,iz);
    ttmp[1] -= d2 * t2m * F3S(p_f, _RMASK, ix,iy,iz);
    F3S(p_f, _RESIS, ix,iy,iz) += fabsf(d1+d2) * F3S(p_f, _ZMASK, ix,iy,iz);
    F3S(p_f, _FLX + XX, ix,iy,iz) = ttmp[0] - ttmp[1];
  } fld3d_foreach_end;
}

static void
bcthy3z_const(fld3d_t p_f, int XX, int YY, int ZZ, int IX, int IY, int IZ,
	      int JX1, int JY1, int JZ1, int JX2, int JY2, int JZ2, mrc_fld_data_t dt, int m_curr)
{
  calc_avg_dz_By(p_f, m_curr, XX, YY, ZZ, JX1, JY1, JZ1, JX2, JY2, JZ2);

  // edge centered E = - v x B (+ dissipation)
  fld3d_foreach(ix,iy,iz, 1, 0) {
    mrc_fld_data_t ttmp[2];
    calc_v_x_B(ttmp, p_f, m_curr, ix, iy, iz, XX, YY, ZZ, IX, IY, IZ,
	       JX1, JY1, JZ1, JX2, JY2, JZ2, dt);

    mrc_fld_data_t vcurrXX = CC_TO_EC(p_f, _CURRX + XX, ix,iy,iz, IX,IY,IZ);
    mrc_fld_data_t vresis = CC_TO_EC(p_f, _RESIS, ix,iy,iz, IX,IY,IZ);
    F3S(p_f, _FLX + XX, ix,iy,iz) = ttmp[0] - ttmp[1] - vresis * vcurrXX;
  } fld3d_foreach_end;
}

static void
calce_nl1_c(fld3d_t p_f, mrc_fld_data_t dt, int m_curr)
{
  bcthy3z_NL1(p_f, 0,1,2, 0,1,1, 0,1,0, 0,0,1, dt, m_curr);
  bcthy3z_NL1(p_f, 1,2,0, 1,0,1, 0,0,1, 1,0,0, dt, m_curr);
  bcthy3z_NL1(p_f, 2,0,1, 1,1,0, 1,0,0, 0,1,0, dt, m_curr);
}

static void
calce_const_c(fld3d_t p_f, mrc_fld_data_t dt, int m_curr)
{
  bcthy3z_const(p_f, 0,1,2, 0,1,1, 0,1,0, 0,0,1, dt, m_curr);
  bcthy3z_const(p_f, 1,2,0, 1,0,1, 0,0,1, 1,0,0, dt, m_curr);
  bcthy3z_const(p_f, 2,0,1, 1,1,0, 1,0,0, 0,1,0, dt, m_curr);
}

// ----------------------------------------------------------------------
// patch_calce_c

static void
patch_calce_c(fld3d_t p_f, mrc_fld_data_t dt, int m_curr)
{
  switch (s_magdiffu) {
  case MAGDIFFU_NL1:
    return calce_nl1_c(p_f, dt, m_curr);
  case MAGDIFFU_CONST:
    return calce_const_c(p_f, dt, m_curr);
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
patch_calce_fortran(fld3d_t p_f, mrc_fld_data_t dt, int m_curr)
{
  calce_F77(F(p_f, _B1X + m_curr), F(p_f, _B1Y + m_curr), F(p_f, _B1Z + m_curr), 
	    F(p_f, _ZMASK), F(p_f, _RMASK), F(p_f, _RESIS), 
	    F(p_f, _FLX), F(p_f, _FLY), F(p_f, _FLZ), 
	    F(p_f, _VX), F(p_f, _VY), F(p_f, _VZ), 
	    F(p_f, _CURRX), F(p_f, _CURRY), F(p_f, _CURRZ), 
	    &dt, &s_mhd_time);
}

#endif

// ----------------------------------------------------------------------
// patch_calce

static void
patch_calce(fld3d_t p_f, mrc_fld_data_t dt, int m_curr)
{
  if (s_opt_mhd_calce == OPT_MHD_C) {
    patch_calce_c(p_f, dt, m_curr);
#if defined(HAVE_OPENGGCM_FORTRAN) && defined(MRC_FLD_AS_FLOAT_H)
  } else if (s_opt_mhd_calce == OPT_MHD_FORTRAN) {
    patch_calce_fortran(p_f, dt, m_curr);
#endif
  } else {
    assert(0);
  }
}

#endif
