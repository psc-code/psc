
// ----------------------------------------------------------------------
// curr_c
//
// edge centered current density

static void
curr_c(fld3d_t p_f, int m_j, int m_curr)
{
  fld3d_foreach(ix,iy,iz, 2, 1) {
    F3S(p_f, m_j + 0, ix,iy,iz) =
      (F3S(p_f, m_curr + _B1Z, ix,iy+1,iz) - F3S(p_f, m_curr + _B1Z, ix,iy,iz)) * BD4Y(iy) -
      (F3S(p_f, m_curr + _B1Y, ix,iy,iz+1) - F3S(p_f, m_curr + _B1Y, ix,iy,iz)) * BD4Z(iz);
    F3S(p_f, m_j + 1, ix,iy,iz) =
      (F3S(p_f, m_curr + _B1X, ix,iy,iz+1) - F3S(p_f, m_curr + _B1X, ix,iy,iz)) * BD4Z(iz) -
      (F3S(p_f, m_curr + _B1Z, ix+1,iy,iz) - F3S(p_f, m_curr + _B1Z, ix,iy,iz)) * BD4X(ix);
    F3S(p_f, m_j + 2, ix,iy,iz) =
      (F3S(p_f, m_curr + _B1Y, ix+1,iy,iz) - F3S(p_f, m_curr + _B1Y, ix,iy,iz)) * BD4X(ix) -
      (F3S(p_f, m_curr + _B1X, ix,iy+1,iz) - F3S(p_f, m_curr + _B1X, ix,iy,iz)) * BD4Y(iy);
  } fld3d_foreach_end;
}

// ----------------------------------------------------------------------
// currbb_c
//
// cell-averaged B

static void
currbb_c(fld3d_t p_f, int m, int m_curr)
{
  fld3d_foreach(ix,iy,iz, 1, 1) {
    F3S(p_f, m+0, ix,iy,iz) = .5f * (F3S(p_f, m_curr + _B1X, ix  ,iy,iz) +
				     F3S(p_f, m_curr + _B1X, ix-1,iy,iz));
    F3S(p_f, m+1, ix,iy,iz) = .5f * (F3S(p_f, m_curr + _B1Y, ix,iy  ,iz) +
				     F3S(p_f, m_curr + _B1Y, ix,iy-1,iz));
    F3S(p_f, m+2, ix,iy,iz) = .5f * (F3S(p_f, m_curr + _B1Z, ix,iy,iz  ) +
				     F3S(p_f, m_curr + _B1Z, ix,iy,iz-1));
  } fld3d_foreach_end;
}

// ----------------------------------------------------------------------
// patch_push_ej_c

static void
patch_push_ej_c(fld3d_t p_f, mrc_fld_data_t dt, int m_curr, int m_next)
{
  enum { XJX = _BX, XJY = _BY, XJZ = _BZ };
  enum { BX = _TMP1, BY = _TMP2, BZ = _TMP3 };

  curr_c(p_f, XJX, m_curr);
  currbb_c(p_f, BX, m_curr);
	
  mrc_fld_data_t s1 = .25f * dt;
  fld3d_foreach(ix,iy,iz, 0, 0) {
    mrc_fld_data_t z = F3S(p_f,_ZMASK, ix,iy,iz);
    mrc_fld_data_t s2 = s1 * z;
    mrc_fld_data_t cx = (F3S(p_f, XJX, ix  ,iy  ,iz  ) +
		F3S(p_f, XJX, ix  ,iy-1,iz  ) +
		F3S(p_f, XJX, ix  ,iy  ,iz-1) +
		F3S(p_f, XJX, ix  ,iy-1,iz-1));
    mrc_fld_data_t cy = (F3S(p_f, XJY, ix  ,iy  ,iz  ) +
		F3S(p_f, XJY, ix-1,iy  ,iz  ) +
		F3S(p_f, XJY, ix  ,iy  ,iz-1) +
		F3S(p_f, XJY, ix-1,iy  ,iz-1));
    mrc_fld_data_t cz = (F3S(p_f, XJZ, ix  ,iy  ,iz  ) +
		F3S(p_f, XJZ, ix-1,iy  ,iz  ) +
		F3S(p_f, XJZ, ix  ,iy-1,iz  ) +
		F3S(p_f, XJZ, ix-1,iy-1,iz  ));
    mrc_fld_data_t ffx = s2 * (cy * F3S(p_f, BZ, ix,iy,iz) -
		      cz * F3S(p_f, BY, ix,iy,iz));
    mrc_fld_data_t ffy = s2 * (cz * F3S(p_f, BX, ix,iy,iz) -
		      cx * F3S(p_f, BZ, ix,iy,iz));
    mrc_fld_data_t ffz = s2 * (cx * F3S(p_f, BY, ix,iy,iz) -
		      cy * F3S(p_f, BX, ix,iy,iz));
    mrc_fld_data_t duu = (ffx * F3S(p_f, _VX, ix,iy,iz) +
		 ffy * F3S(p_f, _VY, ix,iy,iz) +
		 ffz * F3S(p_f, _VZ, ix,iy,iz));

    F3S(p_f, m_next + _RV1X, ix,iy,iz) += ffx;
    F3S(p_f, m_next + _RV1Y, ix,iy,iz) += ffy;
    F3S(p_f, m_next + _RV1Z, ix,iy,iz) += ffz;
    F3S(p_f, m_next + _UU1 , ix,iy,iz) += duu;
  } fld3d_foreach_end;
}

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

// ----------------------------------------------------------------------
// patch_push_ej

static void
patch_push_ej(fld3d_t p_f, mrc_fld_data_t dt, int m_curr, int m_next)
{
  if (s_opt_mhd_push_ej == OPT_MHD_C) {
    patch_push_ej_c(p_f, dt, m_curr, m_next);
#if defined(HAVE_OPENGGCM_FORTRAN) && defined(MRC_FLD_AS_FLOAT_H)
  } else if (s_opt_mhd_push_ej == OPT_MHD_FORTRAN) {
    patch_push_ej_fortran(p_f, dt, m_curr, m_next);
#endif
  } else {
    assert(0);
  }
}
