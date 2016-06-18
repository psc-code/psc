
// ----------------------------------------------------------------------
// patch_pfie3_c

static void
patch_pfie3_c(fld3d_t p_f, mrc_fld_data_t dt,
	      int m_prev, int m_curr, int m_next)
{
  patch_calce(p_f, dt, m_curr);
  patch_bpush1(p_f, dt, m_prev, m_next);
}

// ----------------------------------------------------------------------
// patch_pfie3_fortran

#if defined(HAVE_OPENGGCM_FORTRAN) && defined(MRC_FLD_AS_FLOAT_H)

#define pfie3_F77 F77_FUNC(pfie3,PFIE3)

void pfie3_F77(real *b1x, real *b1y, real *b1z,
	       real *b2x, real *b2y, real *b2z,
	       real *b3x, real *b3y, real *b3z,
	       real *rvx, real *rvy, real *rvz, real *uu,
	       real *zmask, real *rmask, real *resis,
	       real *flx, real *fly, real *flz,
	       real *vx, real *vy, real *vz,
	       real *currx, real *curry, real *currz,
	       real *dt, real *time);

static void
patch_pfie3_fortran(fld3d_t p_f, mrc_fld_data_t dt,
		    int m_prev, int m_curr, int m_next)
{
  pfie3_F77(F(p_f, _B1X + m_prev), F(p_f, _B1Y + m_prev), F(p_f, _B1Z + m_prev),
  	    F(p_f, _B1X + m_curr), F(p_f, _B1Y + m_curr), F(p_f, _B1Z + m_curr),
  	    F(p_f, _B1X + m_next), F(p_f, _B1Y + m_next), F(p_f, _B1Z + m_next),
  	    F(p_f, _RV1X + m_next), F(p_f, _RV1Y + m_next), F(p_f, _RV1Z + m_next),
  	    F(p_f, _UU1 + m_next),
  	    F(p_f, _ZMASK), F(p_f, _RMASK), F(p_f, _RESIS),
  	    F(p_f, _FLX), F(p_f, _FLY), F(p_f, _FLZ),
  	    F(p_f, _VX), F(p_f, _VY), F(p_f, _VZ),
  	    F(p_f, _CURRX), F(p_f, _CURRY), F(p_f, _CURRZ),
  	    &dt, &s_mhd_time);
}

#endif

// ----------------------------------------------------------------------
// patch_pfie3

static void
patch_pfie3(fld3d_t p_f, mrc_fld_data_t dt, int m_prev, int m_curr, int m_next)
{
  if (s_opt_mhd_pfie3 == OPT_MHD_C) {
    patch_pfie3_c(p_f, dt, m_prev, m_curr, m_next);
#if defined(HAVE_OPENGGCM_FORTRAN) && defined(MRC_FLD_AS_FLOAT_H)
  } else if (s_opt_mhd_pfie3 == OPT_MHD_FORTRAN) {
    patch_pfie3_fortran(p_f, dt, m_prev, m_curr, m_next);
#endif
  } else {
    assert(0);
  }
}

