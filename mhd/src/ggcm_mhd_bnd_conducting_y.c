#ifndef __GGCM_MHD_BND_CONDUCTING_Y_C
#define __GGCM_MHD_BND_CONDUCTING_Y_C

#include "ggcm_mhd_bnd_private.h"

#include "ggcm_mhd_defs.h"
#include "ggcm_mhd_private.h"
#include "ggcm_mhd_crds.h"
#include "ggcm_mhd_step.h"

#include <mrc_domain.h>
#include <string.h>
#include <assert.h>

// prototypes for type specific implementations... there must
// be a better way than this
void ggcm_mhd_bnd_conducting_y_fill_ghost_float(struct ggcm_mhd_bnd *bnd,
						struct mrc_fld *fld_base,
						float bntim);
void ggcm_mhd_bnd_conducting_y_fill_ghost_double(struct ggcm_mhd_bnd *bnd,
						 struct mrc_fld *fld_base,
						 float bntim);
void ggcm_mhd_bnd_conducting_y_fill_ghost_gkeyll(struct ggcm_mhd_bnd *bnd,
						 struct mrc_fld *fld_base,
						 float bntim);

struct ggcm_mhd_bnd_conducting_y {
  void (*conducting_y_fill_ghosts)(struct ggcm_mhd_bnd *bnd,
				   struct mrc_fld *fld_base,
				   float bntim);
};
#define ggcm_mhd_bnd_conducting_y(bnd) mrc_to_subobj(bnd, struct ggcm_mhd_bnd_conducting_y)

static void
ggcm_mhd_bnd_conducting_y_setup(struct ggcm_mhd_bnd *bnd)
{
  struct ggcm_mhd_bnd_conducting_y *sub = ggcm_mhd_bnd_conducting_y(bnd);
  
  const char * step_fld_type = mrc_fld_type(bnd->mhd->fld);
  int mhd_type;
  mrc_fld_get_param_int(bnd->mhd->fld, "mhd_type", &mhd_type);
  
  // setup the fld_type dispatch
  if (strcmp(step_fld_type, "float") == 0) {
    sub->conducting_y_fill_ghosts = ggcm_mhd_bnd_conducting_y_fill_ghost_float;
  } else if (strcmp(step_fld_type, "double") == 0) {
    if (mhd_type == MT_GKEYLL) {
      assert(bnd->mhd->fld->_aos);
      assert(bnd->mhd->fld->_c_order);
      sub->conducting_y_fill_ghosts = ggcm_mhd_bnd_conducting_y_fill_ghost_gkeyll;
    }
    else
      sub->conducting_y_fill_ghosts = ggcm_mhd_bnd_conducting_y_fill_ghost_double;
  } else {
    mprintf(">> fld_type = '%s'\n", step_fld_type);
    assert(0); // conducting_y not implemented for this step type
  }
}

static void
_ggcm_mhd_bnd_conducting_y_fill_ghosts(struct ggcm_mhd_bnd *bnd,
				       struct mrc_fld *fld_base,
				       float bntim)
{
  struct ggcm_mhd_bnd_conducting_y *sub = ggcm_mhd_bnd_conducting_y(bnd);
  sub->conducting_y_fill_ghosts(bnd, fld_base, bntim);
}


// ----------------------------------------------------------------------
// ggcm_mhd_bnd_conducting_y2_ops

struct ggcm_mhd_bnd_ops ggcm_mhd_bnd_conducting_y_ops = {
  .name        = "conducting_y",
  .size        = sizeof(struct ggcm_mhd_bnd_conducting_y),
  .setup       = ggcm_mhd_bnd_conducting_y_setup,
  .fill_ghosts = _ggcm_mhd_bnd_conducting_y_fill_ghosts,
};

#endif // __GGCM_MHD_BND_CONDUCTING_Y
