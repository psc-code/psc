
#ifndef PDE_MHD_PRIMBB_C
#define PDE_MHD_PRIMBB_C

#if OPT_STAGGER == OPT_STAGGER_GGCM

#define BTXcc(p_B, i,j,k) (.5f * (F3S(p_B, 0, i,j,k) + F3S(p_B, 0, i-1,j,k)))
#define BTYcc(p_B, i,j,k) (.5f * (F3S(p_B, 1, i,j,k) + F3S(p_B, 1, i,j-1,k)))
#define BTZcc(p_B, i,j,k) (.5f * (F3S(p_B, 2, i,j,k) + F3S(p_B, 2, i,j,k-1)))

#else

#define BTXcc(p_B, i,j,k) (.5f * (BT_(p_B, 0, i,j,k) + BT_(p_B, 0, i+di,j,k)))
#define BTYcc(p_B, i,j,k) (.5f * (BT_(p_B, 1, i,j,k) + BT_(p_B, 1, i,j+dj,k)))
#define BTZcc(p_B, i,j,k) (.5f * (BT_(p_B, 2, i,j,k) + BT_(p_B, 2, i,j,k+dk)))

#endif

// ----------------------------------------------------------------------
// patch_primbb_c
//
// cell-averaged Btotal (ie., add B0 back in, if applicable)
// was also known in Fortran as currbb()

static void
patch_primbb_c(fld3d_t p_Bcc, fld3d_t p_U)
{
  fld3d_t p_B = fld3d_make_view(p_U, BX);

  // FIXME, at least in most cases (zmaskn, get_dt), (0, 0) limits would
  // be enough
  fld3d_foreach(i,j,k, 1, 1) {
    F3S(p_Bcc, 0, i,j,k) = BTXcc(p_B, i,j,k);
    F3S(p_Bcc, 1, i,j,k) = BTYcc(p_B, i,j,k);
    F3S(p_Bcc, 2, i,j,k) = BTZcc(p_B, i,j,k);
  } fld3d_foreach_end;
}

// ----------------------------------------------------------------------
// patch_primbb_fortran

#if defined(HAVE_OPENGGCM_FORTRAN) && defined(MRC_FLD_AS_FLOAT_H)

#include "pde/pde_fortran.h"

#define primbb_F77 F77_FUNC(primbb,PRIMBB)

void primbb_F77(real *bx1, real *by1, real *bz1, real *bx, real *by, real *bz);

static void
patch_primbb_fortran(fld3d_t p_bcc, fld3d_t p_U)
{
  primbb_F77(F(p_U, BX), F(p_U, BY), F(p_U, BZ),
	     F(p_bcc, 0), F(p_bcc, 1), F(p_bcc, 2));
}

#endif

// ----------------------------------------------------------------------
// patch_primbb

static void _mrc_unused
patch_primbb(fld3d_t p_bcc, fld3d_t p_U)
{
  if (s_opt_mhd_primbb == OPT_MHD_C) {
    patch_primbb_c(p_bcc, p_U);
#if defined(HAVE_OPENGGCM_FORTRAN) && defined(MRC_FLD_AS_FLOAT_H)
  } else if (s_opt_mhd_primbb == OPT_MHD_FORTRAN) {
    patch_primbb_fortran(p_bcc, p_U);
#endif
  } else {
    assert(0);
  }
}

#endif
