
// ----------------------------------------------------------------------
// patch_zmaskn_c

static void
patch_zmaskn_c(fld3d_t p_f)
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

// ----------------------------------------------------------------------
// patch_zmaskn_fortran

#if defined(HAVE_OPENGGCM_FORTRAN) && defined(MRC_FLD_AS_FLOAT_H)

#include "pde/pde_fortran.h"

#define zmaskn_F77 F77_FUNC(zmaskn,ZMASKN)

void zmaskn_F77(real *rr, real *bx, real *by, real *bz, 
		real *zmask, real *ymask);

static void
patch_zmaskn_fortran(fld3d_t p_f)
{
  zmaskn_F77(F(p_f, _RR), F(p_f, _BX), F(p_f, _BY), F(p_f, _BZ),
	     F(p_f, _ZMASK), F(p_f, _YMASK));
}

#endif

// ----------------------------------------------------------------------
// patch_zmaskn

static void _mrc_unused
patch_zmaskn(fld3d_t p_f)
{
  if (s_opt_mhd_zmaskn == OPT_MHD_C) {
    patch_zmaskn_c(p_f);
#if defined(HAVE_OPENGGCM_FORTRAN) && defined(MRC_FLD_AS_FLOAT_H)
  } else if (s_opt_mhd_zmaskn == OPT_MHD_FORTRAN) {
    patch_zmaskn_fortran(p_f);
#endif
  } else {
    assert(0);
  }
}

