
#include "ggcm_mhd_defs.h"
#include "ggcm_mhd_defs_extra.h"
#include "ggcm_mhd_private.h"

#include <mrc_profile.h>

#include <math.h>
#include <assert.h>

// FIXME: major ugliness
// The fortran fields do primitive vars in the order _RR,_PP,_VX,_VY,_VZ
// but in C, we stick with the corresponding conservative var order, ie.,
// RR,VX,VY,VZ,PP
// The below hackily switches the order around in C, so that it matches fortran

#define PP 1
#define VX 2
#define VY 3
#define VZ 4

#include "pde/pde_defs.h"

// mhd options

#define OPT_EQN OPT_EQN_MHD_SCONS

#include "pde/pde_setup.c"
#include "pde/pde_mhd_setup.c"
#include "pde/pde_mhd_zmaskn.c"

// ----------------------------------------------------------------------
// zmaskn_c

void
zmaskn_c(struct ggcm_mhd *mhd)
{
  static int PR;
  if (!PR) {
    PR = prof_register("zmaskn_c", 1., 0, 0);
  }
  prof_start(PR);

  struct mrc_fld *f = mrc_fld_get_as(mhd->fld, FLD_TYPE);

  pde_setup(f);
  pde_mhd_setup(mhd);

  fld3d_t p_f;
  fld3d_t p_zmask, p_W, p_bcc, p_ymask;
  fld3d_setup(&p_f, f);

  for (int p = 0; p < mrc_fld_nr_patches(f); p++) {
    fld3d_get(&p_f, p);
    fld3d_setup_view(&p_zmask, p_f, _ZMASK);
    fld3d_setup_view(&p_W    , p_f, _RR);
    fld3d_setup_view(&p_bcc  , p_f, _BX);
    fld3d_setup_view(&p_ymask, p_f, _YMASK);

    patch_zmaskn(p_zmask, p_W, p_bcc, p_ymask);

    fld3d_put(&p_f, p);
  }
  mrc_fld_put_as(f, mhd->fld);
  prof_stop(PR);
}

