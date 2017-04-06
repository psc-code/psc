
#ifndef PDE_MHD_CALC_CURRENT_C
#define PDE_MHD_CALC_CURRENT_C

#include "pde/pde_mhd_setup.c"

// ----------------------------------------------------------------------
// patch_calc_current_ec
//
// (original Fortran name: curr())
// edge centered current density

#if OPT_STAGGER == OPT_STAGGER_GGCM

static void _mrc_unused
patch_calc_current_ec(fld3d_t p_J, fld3d_t p_B)
{
  fld3d_foreach(i,j,k, 2, 1) {
    F3S(p_J, 0, i,j,k) = s_mu0_inv * ((F3S(p_B, 2, i,j+dj,k) - F3S(p_B, 2, i,j,k)) * PDE_INV_DYF(j+1) -
				      (F3S(p_B, 1, i,j,k+dk) - F3S(p_B, 1, i,j,k)) * PDE_INV_DZF(k+1));
    F3S(p_J, 1, i,j,k) = s_mu0_inv * ((F3S(p_B, 0, i,j,k+dk) - F3S(p_B, 0, i,j,k)) * PDE_INV_DZF(k+1) -
				      (F3S(p_B, 2, i+di,j,k) - F3S(p_B, 2, i,j,k)) * PDE_INV_DXF(i+1));
    F3S(p_J, 2, i,j,k) = s_mu0_inv * ((F3S(p_B, 1, i+di,j,k) - F3S(p_B, 1, i,j,k)) * PDE_INV_DXF(i+1) -
				      (F3S(p_B, 0, i,j+dj,k) - F3S(p_B, 0, i,j,k)) * PDE_INV_DYF(j+1));
  } fld3d_foreach_end;
}

#else

static void
patch_calc_current_ec(fld3d_t p_J, fld3d_t p_B)
{
  fld3d_foreach(i,j,k, 1, 2) {
    F3S(p_J, 0, i,j,k) = s_mu0_inv * ((F3S(p_B, 2, i,j,k) - F3S(p_B, 2, i,j-dj,k)) * PDE_INV_DYF(j) -
				      (F3S(p_B, 1, i,j,k) - F3S(p_B, 1, i,j,k-dk)) * PDE_INV_DZF(k));
    F3S(p_J, 1, i,j,k) = s_mu0_inv * ((F3S(p_B, 0, i,j,k) - F3S(p_B, 0, i,j,k-dk)) * PDE_INV_DZF(k) -
				      (F3S(p_B, 2, i,j,k) - F3S(p_B, 2, i-di,j,k)) * PDE_INV_DXF(i));
    F3S(p_J, 2, i,j,k) = s_mu0_inv * ((F3S(p_B, 1, i,j,k) - F3S(p_B, 1, i-di,j,k)) * PDE_INV_DXF(i) -
				      (F3S(p_B, 0, i,j,k) - F3S(p_B, 0, i,j-dj,k)) * PDE_INV_DYF(j));
  } fld3d_foreach_end;
}

#endif

// ----------------------------------------------------------------------
// patch_calc_current_cc
//
// (original Fortran name: curbc())
// cell centered current density

static void _mrc_unused
patch_calc_current_cc(fld3d_t p_Jcc, fld3d_t p_U, fld3d_t p_zmask)
{ 
  static fld3d_t p_Jec;
  fld3d_setup_tmp_compat(&p_Jec, 3, _TMP1); /* was named _TX */
  fld3d_t p_B = fld3d_make_view(p_U, BX);

  patch_calc_current_ec(p_Jec, p_B);

  // j averaged to cell-centered
  // FIXME, multiplication by z_mask is redundant?
  fld3d_foreach(i,j,k, 1, 1) {
    mrc_fld_data_t z = F3S(p_zmask, 0, i, j, k);
    F3S(p_Jcc, 0, i,j,k) = z * EC_TO_CC(p_Jec, 0, i,j,k);
    F3S(p_Jcc, 1, i,j,k) = z * EC_TO_CC(p_Jec, 1, i,j,k);
    F3S(p_Jcc, 2, i,j,k) = z * EC_TO_CC(p_Jec, 2, i,j,k);
  } fld3d_foreach_end;
}

// ----------------------------------------------------------------------
// patch_calc_current_cc_bgrid_fc
//
// doesn't include zmask

static void _mrc_unused
patch_calc_current_cc_bgrid_fc(fld3d_t p_Jcc, fld3d_t p_B)
{ 
  static fld3d_t p_Jec;
  fld3d_setup_tmp(&p_Jec, 3);

  patch_calc_current_ec(p_Jec, p_B);

  // j averaged to cell-centered
  fld3d_foreach(i,j,k, 1, 1) {
    F3S(p_Jcc, 0, i,j,k) = EC_TO_CC(p_Jec, 0, i,j,k);
    F3S(p_Jcc, 1, i,j,k) = EC_TO_CC(p_Jec, 1, i,j,k);
    F3S(p_Jcc, 2, i,j,k) = EC_TO_CC(p_Jec, 2, i,j,k);
  } fld3d_foreach_end;
}

// ----------------------------------------------------------------------
// patch_calc_current_cc_bgrid_cc
//
// calculate cell-centered J from cell-centered B

static void _mrc_unused
patch_calc_current_cc_bgrid_cc(fld3d_t p_J, fld3d_t p_B)
{
  fld3d_foreach(i,j,k, 1, 1) {
    F3S(p_J, 0, i,j,k) = s_mu0_inv *
      ((F3S(p_B, 2, i,j+dj,k) - F3S(p_B, 2, i,j-dj,k)) * .5f * PDE_INV_DY(j) -
       (F3S(p_B, 1, i,j,k+dk) - F3S(p_B, 1, i,j,k-dk)) * .5f * PDE_INV_DZ(k));
    F3S(p_J, 1, i,j,k) = s_mu0_inv *
      ((F3S(p_B, 0, i,j,k+dk) - F3S(p_B, 0, i,j,k-dk)) * .5f * PDE_INV_DZ(k) -
       (F3S(p_B, 2, i+di,j,k) - F3S(p_B, 2, i-di,j,k)) * .5f * PDE_INV_DX(i));
    F3S(p_J, 2, i,j,k) = s_mu0_inv *
      ((F3S(p_B, 1, i+di,j,k) - F3S(p_B, 1, i-di,j,k)) * .5f * PDE_INV_DX(i) -
       (F3S(p_B, 0, i,j+dj,k) - F3S(p_B, 0, i,j-dj,k)) * .5f * PDE_INV_DY(j));
  } fld3d_foreach_end;
}

#endif
