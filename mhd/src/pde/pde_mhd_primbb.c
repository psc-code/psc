
// ----------------------------------------------------------------------
// patch_primbb

static void
patch_primbb(fld3d_t p_f, int m)
{
  fld3d_foreach(i,j,k, 1, 2) {
    F3S(p_f,_BX, i,j,k) = .5f * (F3S(p_f, m + _B1X, i,j,k) +
				 F3S(p_f, m + _B1X, i-1,j,k));
    F3S(p_f,_BY, i,j,k) = .5f * (F3S(p_f, m + _B1Y, i,j,k) +
				 F3S(p_f, m + _B1Y, i,j-1,k));
    F3S(p_f,_BZ, i,j,k) = .5f * (F3S(p_f, m + _B1Z, i,j,k) +
				 F3S(p_f, m + _B1Z, i,j,k-1));
  } fld3d_foreach_end;
}
