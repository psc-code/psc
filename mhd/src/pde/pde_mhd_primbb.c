
#ifndef PDE_MHD_PRIMBB_C
#define PDE_MHD_PRIMBB_C

// ----------------------------------------------------------------------
// patch_primbb_c
//
// was also known in Fortran as currbb()

#if OPT_STAGGER == OPT_STAGGER_GGCM

#define BXcc(p_U, i,j,k) (.5f*(F3S(p_U, BX, i,j,k) + F3S(p_U, BX,i-1,j,k)))
#define BYcc(p_U, i,j,k) (.5f*(F3S(p_U, BY, i,j,k) + F3S(p_U, BY,i,j-1,k)))
#define BZcc(p_U, i,j,k) (.5f*(F3S(p_U, BZ, i,j,k) + F3S(p_U, BZ,i,j,k-1)))

#else

#define BXcc(p_U, i,j,k) (.5f*(F3S(p_U, BX, i,j,k) + F3S(p_U, BX,i+1,j,k)))
#define BYcc(p_U, i,j,k) (.5f*(F3S(p_U, BY, i,j,k) + F3S(p_U, BY,i,j+1,k)))
#define BZcc(p_U, i,j,k) (.5f*(F3S(p_U, BZ, i,j,k) + F3S(p_U, BZ,i,j,k+1)))

#endif

static void
patch_primbb_c(fld3d_t p_bcc, fld3d_t p_U)
{
  fld3d_foreach(i,j,k, 1, 2) {
    F3S(p_bcc, 0, i,j,k) = BXcc(p_U, i,j,k);
    F3S(p_bcc, 1, i,j,k) = BYcc(p_U, i,j,k);
    F3S(p_bcc, 2, i,j,k) = BZcc(p_U, i,j,k);
  } fld3d_foreach_end;
}

// ----------------------------------------------------------------------
// patch_calc_Bt_cc
//
// cell-averaged Btotal (ie., add B0 back in, if applicable)
// FIXME, consolidate with above

static void _mrc_unused
patch_calc_Bt_cc(fld3d_t p_Bcc, fld3d_t p_U, int l, int r)
{
  fld3d_t p_B = fld3d_make_view(p_U, BX);

  fld3d_foreach(i,j,k, l, r) {
    F3S(p_Bcc, 0, i,j,k) = .5f * (BT_(p_B, 0, i,j,k) + BT_(p_B, 0, i+di,j,k));
    F3S(p_Bcc, 1, i,j,k) = .5f * (BT_(p_B, 1, i,j,k) + BT_(p_B, 1, i,j+dj,k));
    F3S(p_Bcc, 2, i,j,k) = .5f * (BT_(p_B, 2, i,j,k) + BT_(p_B, 2, i,j,k+dk));
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
