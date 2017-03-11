
#ifndef PDE_MHD_DIVB_C
#define PDE_MHD_DIVB_C

#include "pde/pde_mhd_setup.c"

// ----------------------------------------------------------------------
// patch_calc_divb_grid_cc

static void _mrc_unused
patch_calc_divb_bgrid_cc(fld3d_t p_divB, fld3d_t p_B)
{
  fld3d_foreach(i,j,k, 0, 0) {
    F3S(p_divB, 0, i,j,k) = .5f * ((F3S(p_B, 0, i+di,j,k) - F3S(p_B, 0, i-di,j,k)) * PDE_INV_DX(i) +
				   (F3S(p_B, 1, i,j+dj,k) - F3S(p_B, 1, i,j-dj,k)) * PDE_INV_DY(j) +
				   (F3S(p_B, 2, i,j,k+dk) - F3S(p_B, 2, i,j,k-dk)) * PDE_INV_DZ(k));
  } fld3d_foreach_end;
}

// ----------------------------------------------------------------------
// patch_calc_divb_grid_fc

static void _mrc_unused
patch_calc_divb_bgrid_fc(fld3d_t p_divB, fld3d_t p_B)
{
  fld3d_foreach(i,j,k, 0, 0) {
#if OPT_STAGGER == OPT_STAGGER_GGCM
    F3S(p_divB, 0, i,j,k) = ((F3S(p_B, 0, i,j,k) - F3S(p_B, 0, i-di,j,k)) * PDE_INV_DX(i) +
			     (F3S(p_B, 1, i,j,k) - F3S(p_B, 1, i,j-dj,k)) * PDE_INV_DY(j) +
			     (F3S(p_B, 2, i,j,k) - F3S(p_B, 2, i,j,k-dk)) * PDE_INV_DZ(k));
#else
    F3S(p_divB, 0, i,j,k) = ((F3S(p_B, 0, i+di,j,k) - F3S(p_B, 0, i,j,k)) * PDE_INV_DX(i) +
			     (F3S(p_B, 1, i,j+dj,k) - F3S(p_B, 1, i,j,k)) * PDE_INV_DY(j) +
			     (F3S(p_B, 2, i,j,k+dk) - F3S(p_B, 2, i,j,k)) * PDE_INV_DZ(k));
#endif
  } fld3d_foreach_end;
}



#endif
