
#include "pde/pde_defs.h"

#if MT == MT_GRID_FC_GGCM
#define OPT_STAGGER OPT_STAGGER_GGCM
#else
#define OPT_STAGGER OPT_STAGGER_REG
#endif

#include "pde/pde_mhd_calc_current.c"
#include "pde/pde_mhd_divb.c"

// ----------------------------------------------------------------------
// ggcm_mhd_calc_curcc_bgrid_*

void
BGRID_SFX(ggcm_mhd_calc_currcc)(struct ggcm_mhd *mhd, struct mrc_fld *f, int m,
				struct mrc_fld *c)
{
  static bool is_setup = false;
  if (!is_setup) {
    pde_mhd_setup(mhd, mrc_fld_nr_comps(f));
    is_setup = true;
  }
  
  fld3d_t p_U, p_J;
  fld3d_setup(&p_U, f);
  fld3d_setup(&p_J, c);
    
  for (int p = 0; p < mrc_fld_nr_patches(f); p++) {
    pde_patch_set(p);
    fld3d_get(&p_U, p);
    fld3d_get(&p_J, p);
    fld3d_t p_B = fld3d_make_view(p_U, m);

#if MT == MT_BGRID_CC
    patch_calc_current_cc_bgrid_cc(p_J, p_B);
#elif MT == MT_BGRID_FC || MT == MT_BGRID_FC_GGCM
    patch_calc_current_cc_bgrid_fc(p_J, p_B);
#else
#error unknown MT
#endif
  }

  if (0) {
    pde_free();
  }
}

void
BGRID_SFX(ggcm_mhd_calc_divb)(struct ggcm_mhd *mhd, struct mrc_fld *f, struct mrc_fld *divB)
{
  static bool is_setup = false;
  if (!is_setup) {
    pde_mhd_setup(mhd, mrc_fld_nr_comps(f));
    is_setup = true;
  }
  
  fld3d_t p_U, p_divB;
  fld3d_setup(&p_U, f);
  fld3d_setup(&p_divB, divB);

  for (int p = 0; p < mrc_fld_nr_patches(f); p++) {
    pde_patch_set(p);
    fld3d_get(&p_U, p);
    fld3d_get(&p_divB, p);

#if MT == MT_BGRID_CC
    patch_calc_divb_bgrid_cc(p_divB, fld3d_make_view(p_U, BX));
#elif MT == MT_BGRID_FC || MT == MT_BGRID_FC_GGCM
    patch_calc_divb_bgrid_fc(p_divB, fld3d_make_view(p_U, BX));
#else
#error unknown MT
#endif
  }

  if (0) {
    pde_free();
  }
}

