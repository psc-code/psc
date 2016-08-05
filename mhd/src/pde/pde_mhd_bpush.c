
#ifndef PDE_MHD_BPUSH_C
#define PDE_MHD_BPUSH_C

// ----------------------------------------------------------------------
// patch_bpush1_c

static void
patch_bpush1_c(fld3d_t p_f, mrc_fld_data_t dt, int m_prev, int m_next)
{
  fld3d_foreach(ix,iy,iz, 0, 0) {
    F3S(p_f, m_next + _B1X, ix,iy,iz) = F3S(p_f, m_prev + _B1X, ix,iy,iz) +
      dt * (BD3Y(iy) * (F3S(p_f,_FLZ, ix,iy,iz) - F3S(p_f,_FLZ, ix,iy-1,iz)) -
	    BD3Z(iz) * (F3S(p_f,_FLY, ix,iy,iz) - F3S(p_f,_FLY, ix,iy,iz-1)));
    F3S(p_f, m_next + _B1Y, ix,iy,iz) = F3S(p_f, m_prev + _B1Y, ix,iy,iz) +
      dt * (BD3Z(iz) * (F3S(p_f,_FLX, ix,iy,iz) - F3S(p_f,_FLX, ix,iy,iz-1)) -
	    BD3X(ix) * (F3S(p_f,_FLZ, ix,iy,iz) - F3S(p_f,_FLZ, ix-1,iy,iz)));
    F3S(p_f, m_next + _B1Z, ix,iy,iz) = F3S(p_f, m_prev + _B1Z, ix,iy,iz) +
      dt * (BD3X(ix) * (F3S(p_f,_FLY, ix,iy,iz) - F3S(p_f,_FLY, ix-1,iy,iz)) -
	    BD3Y(iy) * (F3S(p_f,_FLX, ix,iy,iz) - F3S(p_f,_FLX, ix,iy-1,iz)));
  } fld3d_foreach_end;
}

// ----------------------------------------------------------------------
// patch_bpush1_fortran

#if defined(HAVE_OPENGGCM_FORTRAN) && defined(MRC_FLD_AS_FLOAT_H)

#define bpush1_F77 F77_FUNC(bpush1,BPUSH1)

void bpush1_F77(real *bx1, real *by1, real *bz1,
		real *bx3, real *by3, real *bz3,
		real *flx, real *fly, real *flz, real *dt);

static void
patch_bpush1_fortran(fld3d_t p_f, mrc_fld_data_t dt, int m_prev, int m_next)
{
  bpush1_F77(F(p_f, _B1X + m_prev), F(p_f, _B1Y + m_prev), F(p_f, _B1Z + m_prev), 
	     F(p_f, _B1X + m_next), F(p_f, _B1Y + m_next), F(p_f, _B1Z + m_next), 
	     F(p_f, _FLX), F(p_f, _FLY), F(p_f, _FLZ), 
	     &dt);
}

#endif

// ----------------------------------------------------------------------
// patch_bpush1

static void
patch_bpush1(fld3d_t p_f, mrc_fld_data_t dt, int m_prev, int m_next)
{
  if (s_opt_mhd_bpush1 == OPT_MHD_C) {
    patch_bpush1_c(p_f, dt, m_prev, m_next);
#if defined(HAVE_OPENGGCM_FORTRAN) && defined(MRC_FLD_AS_FLOAT_H)
  } else if (s_opt_mhd_bpush1 == OPT_MHD_FORTRAN) {
    patch_bpush1_fortran(p_f, dt, m_prev, m_next);
#endif
  } else {
    assert(0);
  }
}

#endif
