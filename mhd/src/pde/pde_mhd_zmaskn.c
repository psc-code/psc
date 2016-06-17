
// ----------------------------------------------------------------------
// patch_zmaskn

static void
patch_zmaskn(fld3d_t p_f)
{
  mrc_fld_data_t va02i = 1.f / sqr(s_speedlimit_code);
  mrc_fld_data_t eps = 1e-15f;

  fld3d_foreach(ix,iy,iz, 2, 2) {
    float bb = 
      sqr(F3S(p_f,_BX, ix,iy,iz)) + 
      sqr(F3S(p_f,_BY, ix,iy,iz)) +
      sqr(F3S(p_f,_BZ, ix,iy,iz));
    float rrm = mrc_fld_max(eps, bb * va02i);
    F3S(p_f, _ZMASK, ix,iy,iz) = F3S(p_f, _YMASK, ix,iy,iz) * 
      mrc_fld_min(1.f, F3S(p_f, _RR, ix,iy,iz) / rrm);
  } fld3d_foreach_end;
}

