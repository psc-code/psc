
// ----------------------------------------------------------------------
// patch_curr
//
// edge centered current density

static void _mrc_unused
patch_curr(fld3d_t p_f, int m_j, int m_curr)
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
// patch_curbc
//
// cell centered current density

static void _mrc_unused
patch_curbc(fld3d_t p_f, int m_curr)
{ 
  enum { _TX = _TMP1, _TY = _TMP2, _TZ = _TMP3 };

  patch_curr(p_f, _TX, m_curr);

  // j averaged to cell-centered
  fld3d_foreach(ix,iy,iz, 1, 1) {
    mrc_fld_data_t s = .25f * F3S(p_f, _ZMASK, ix, iy, iz);
    F3S(p_f, _CURRX, ix,iy,iz) = s * (F3S(p_f, _TX, ix,iy  ,iz  ) + F3S(p_f, _TX, ix,iy-1,iz  ) +
				      F3S(p_f, _TX, ix,iy  ,iz-1) + F3S(p_f, _TX, ix,iy-1,iz-1));
    F3S(p_f, _CURRY, ix,iy,iz) = s * (F3S(p_f, _TY, ix  ,iy,iz  ) + F3S(p_f, _TY, ix-1,iy,iz  ) +
				      F3S(p_f, _TY, ix  ,iy,iz-1) + F3S(p_f, _TY, ix-1,iy,iz-1));
    F3S(p_f, _CURRZ, ix,iy,iz) = s * (F3S(p_f, _TZ, ix  ,iy  ,iz) + F3S(p_f, _TZ, ix-1,iy  ,iz) +
				      F3S(p_f, _TZ, ix  ,iy-1,iz) + F3S(p_f, _TZ, ix-1,iy-1,iz));
  } fld3d_foreach_end;
}

