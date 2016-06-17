
// ----------------------------------------------------------------------
// patch_primbb_c

static void
patch_primbb_c(fld3d_t p_f, int m)
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

// ----------------------------------------------------------------------
// patch_primbb_fortran

#if defined(HAVE_OPENGGCM_FORTRAN) && defined(MRC_FLD_AS_FLOAT_H)

#include "pde/pde_fortran.h"

#define primbb_F77 F77_FUNC(primbb,PRIMBB)

void primbb_F77(real *bx1, real *by1, real *bz1, real *bx, real *by, real *bz);

static void
patch_primbb_fortran(fld3d_t p_f, int m)
{
  if (m == _RR1) {
    primbb_F77(F(p_f, _B1X), F(p_f, _B1Y), F(p_f, _B1Z),
	       F(p_f, _BX), F(p_f, _BY), F(p_f, _BZ));
  } else if (m == _RR2) {
    primbb_F77(F(p_f, _B2X), F(p_f, _B2Y), F(p_f, _B2Z),
	       F(p_f, _BX), F(p_f, _BY), F(p_f, _BZ));
  } else {
    assert(0);
  }
}

#endif

// ----------------------------------------------------------------------
// patch_primbb

static void _mrc_unused
patch_primbb(fld3d_t p_f, int m)
{
  if (s_opt_mhd_primbb == OPT_MHD_C) {
    patch_primbb_c(p_f, m);
#if defined(HAVE_OPENGGCM_FORTRAN) && defined(MRC_FLD_AS_FLOAT_H)
  } else if (s_opt_mhd_primbb == OPT_MHD_FORTRAN) {
    patch_primbb_fortran(p_f, m);
#endif
  } else {
    assert(0);
  }
}

