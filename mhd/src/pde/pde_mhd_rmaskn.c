
#ifndef PDE_MHD_RMASKN
#define PDE_MHD_RMASKN

// ----------------------------------------------------------------------
// patch_rmaskn_c

static void
patch_rmaskn_c(fld3d_t p_rmask, fld3d_t p_zmask)
{
  mrc_fld_data_t diffco = s_diffco;

  fld3d_foreach(i,j,k, 2, 2) {
    F3S(p_rmask, 0, i,j,k) = 0.f;
    if (PDE_CRDX_CC(i) < s_diff_swbnd)
      continue;
    if (j + s_patch.off[1] < s_diff_obnd)
      continue;
    if (k + s_patch.off[2] < s_diff_obnd)
      continue;
    if (i + s_patch.off[0] >= s_gdims[0] - s_diff_obnd)
      continue;
    if (j + s_patch.off[1] >= s_gdims[1] - s_diff_obnd)
      continue;
    if (k + s_patch.off[2] >= s_gdims[2] - s_diff_obnd)
      continue;
    F3S(p_rmask, 0, i,j,k) = diffco * F3S(p_zmask, 0, i,j,k);
  } fld3d_foreach_end;
}

// ----------------------------------------------------------------------
// patch_rmaskn_fortran

#if defined(HAVE_OPENGGCM_FORTRAN) && defined(MRC_FLD_AS_FLOAT_H)

#include "pde/pde_fortran.h"

#define rmaskn_F77 F77_FUNC(rmaskn,RMASKN)

void rmaskn_F77(real *zmask, real *rmask);

static void
patch_rmaskn_fortran(fld3d_t p_rmask, fld3d_t p_zmask)
{
  rmaskn_F77(F(p_zmask, 0), F(p_rmask, 0));
}

#endif

// ----------------------------------------------------------------------
// patch_rmaskn

static void _mrc_unused
patch_rmaskn(fld3d_t p_rmask, fld3d_t p_zmask)
{
  if (s_opt_mhd_rmaskn == OPT_MHD_C) {
    patch_rmaskn_c(p_rmask, p_zmask);
#if defined(HAVE_OPENGGCM_FORTRAN) && defined(MRC_FLD_AS_FLOAT_H)
  } else if (s_opt_mhd_rmaskn == OPT_MHD_FORTRAN) {
    patch_rmaskn_fortran(p_rmask, p_zmask);
#endif
  } else {
    assert(0);
  }
}

#endif
