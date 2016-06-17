
// ----------------------------------------------------------------------
// patch_rmaskn_c

static void
patch_rmaskn_c(fld3d_t p_f)
{
  mrc_fld_data_t diffco = s_diffco;

  fld3d_foreach(ix,iy,iz, 2, 2) {
    F3S(p_f,_RMASK, ix,iy,iz) = 0.f;
    if (FX1X(ix) < s_diff_swbnd)
      continue;
    if (iy + s_patch.off[1] < s_diff_obnd)
      continue;
    if (iz + s_patch.off[2] < s_diff_obnd)
      continue;
    if (ix + s_patch.off[0] >= s_gdims[0] - s_diff_obnd)
      continue;
    if (iy + s_patch.off[1] >= s_gdims[1] - s_diff_obnd)
      continue;
    if (iz + s_patch.off[2] >= s_gdims[2] - s_diff_obnd)
      continue;
    F3S(p_f, _RMASK, ix,iy,iz) = diffco * F3S(p_f, _ZMASK, ix,iy,iz);
  } fld3d_foreach_end;
}

// ----------------------------------------------------------------------
// patch_rmaskn_fortran

#if defined(HAVE_OPENGGCM_FORTRAN) && defined(MRC_FLD_AS_FLOAT_H)

#include "pde/pde_fortran.h"

#define rmaskn_F77 F77_FUNC(rmaskn,RMASKN)

void rmaskn_F77(real *zmask, real *rmask);

static void
patch_rmaskn_fortran(fld3d_t p_f)
{
  rmaskn_F77(F(p_f, _ZMASK), F(p_f, _RMASK));
}

#endif

// ----------------------------------------------------------------------
// patch_rmaskn

static void _mrc_unused
patch_rmaskn(fld3d_t p_f)
{
  if (s_opt_mhd_rmaskn == OPT_MHD_C) {
    patch_rmaskn_c(p_f);
#if defined(HAVE_OPENGGCM_FORTRAN) && defined(MRC_FLD_AS_FLOAT_H)
  } else if (s_opt_mhd_rmaskn == OPT_MHD_FORTRAN) {
    patch_rmaskn_fortran(p_f);
#endif
  } else {
    assert(0);
  }
}

