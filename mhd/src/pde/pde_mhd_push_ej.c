
// ----------------------------------------------------------------------
// patch_push_ej_fortran

#if defined(HAVE_OPENGGCM_FORTRAN) && defined(MRC_FLD_AS_FLOAT_H)

#define push_ej_F77 F77_FUNC(push_ej,PUSH_EJ)

void push_ej_F77(real *b1x, real *b1y, real *b1z,
		 real *rv1x, real *rv1y, real *rv1z, real *uu1,
		 real *zmask, real *vx, real *vy, real *vz,
		 real *dt);

static void
patch_push_ej_fortran(fld3d_t p_f, mrc_fld_data_t dt, int m_curr, int m_next)
{
  push_ej_F77(F(p_f, _B1X + m_curr), F(p_f, _B1Y + m_curr), F(p_f, _B1Z + m_curr),
	      F(p_f, _RV1X + m_next), F(p_f, _RV1Y + m_next), F(p_f, _RV1Z + m_next),
	      F(p_f, _UU1 + m_next), 
	      F(p_f, _ZMASK), F(p_f, _VX), F(p_f, _VY), F(p_f, _VZ), &dt);
}

#endif
