
#include "ggcm_mhd_private.h"
#include "ggcm_mhd_defs.h"

#include <mrc_fld_as_double.h>

#include "pde/pde_defs.h"

#define OPT_STAGGER OPT_STAGGER_REG
#define OPT_FLD1D OPT_FLD1D_C_ARRAY

#include "pde/pde_mhd_calc_current.c"

// ----------------------------------------------------------------------
// ggcm_mhd_calc_curcc_bgrid_fc

void
ggcm_mhd_calc_currcc_bgrid_fc(struct ggcm_mhd *mhd, struct mrc_fld *f, int m,
			      struct mrc_fld *c)
{
  static bool is_setup = false;
  if (!is_setup) {
    pde_setup(f);
    pde_mhd_setup(mhd);
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

    patch_calc_current_cc_bgrid_fc(p_J, p_B);
  }

  if (0) {
    pde_free();
  }
}

