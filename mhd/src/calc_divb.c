
#include "ggcm_mhd_defs.h"
#include "ggcm_mhd_private.h"
#include "ggcm_mhd_crds.h"

#include <mrc_fld.h>
#include <mrc_domain.h>

#include <stdio.h>
#include <math.h>
#include <assert.h>

#include <mrc_fld_as_double.h>

#include "pde/pde_defs.h"

#define OPT_FLD1D OPT_FLD1D_C_ARRAY

#include "pde/pde_mhd_divb.c"

static void
ggcm_mhd_calc_divb_bgrid_cc(struct ggcm_mhd *mhd, struct mrc_fld *f, struct mrc_fld *divB)
{
  fld3d_t p_U, p_divB;
  fld3d_setup(&p_U, f);
  fld3d_setup(&p_divB, divB);

  for (int p = 0; p < mrc_fld_nr_patches(f); p++) {
    pde_patch_set(p);
    fld3d_get(&p_U, p);
    fld3d_get(&p_divB, p);

    patch_calc_divb_bgrid_cc(p_divB, fld3d_make_view(p_U, BX));
  }
}

static void
ggcm_mhd_calc_divb_bgrid_fc(struct ggcm_mhd *mhd, struct mrc_fld *f, struct mrc_fld *divB)
{
  fld3d_t p_U, p_divB;
  fld3d_setup(&p_U, f);
  fld3d_setup(&p_divB, divB);

  for (int p = 0; p < mrc_fld_nr_patches(f); p++) {
    pde_patch_set(p);
    fld3d_get(&p_U, p);
    fld3d_get(&p_divB, p);

    patch_calc_divb_bgrid_fc(p_divB, fld3d_make_view(p_U, BX));
  }
}

static void
ggcm_mhd_calc_divb_bgrid_fc_ggcm(struct ggcm_mhd *mhd, struct mrc_fld *f, struct mrc_fld *divB)
{
  fld3d_t p_U, p_divB;
  fld3d_setup(&p_U, f);
  fld3d_setup(&p_divB, divB);

  for (int p = 0; p < mrc_fld_nr_patches(f); p++) {
    pde_patch_set(p);
    fld3d_get(&p_U, p);
    fld3d_get(&p_divB, p);

    patch_calc_divb_bgrid_fc_ggcm(p_divB, fld3d_make_view(p_U, BX));
  }
}

// ----------------------------------------------------------------------
// ggcm_mhd_calc_divb

void
ggcm_mhd_calc_divb(struct ggcm_mhd *mhd, struct mrc_fld *fld, struct mrc_fld *divb)
{
  static bool is_setup = false;
  if (!is_setup) {
    pde_setup(fld);
    pde_mhd_setup(mhd);
    is_setup = true;
  }

  int mhd_type;
  mrc_fld_get_param_int(fld, "mhd_type", &mhd_type);

  struct mrc_fld *f = mrc_fld_get_as(fld, FLD_TYPE);
  struct mrc_fld *d = mrc_fld_get_as(divb, FLD_TYPE);

  if (MT_BGRID(mhd_type) == MT_BGRID_CC) {
    ggcm_mhd_calc_divb_bgrid_cc(mhd, f, d);
  } else if (MT_BGRID(mhd_type) == MT_BGRID_FC) {
    ggcm_mhd_calc_divb_bgrid_fc(mhd, f, d);
  } else if (MT_BGRID(mhd_type) == MT_BGRID_FC_GGCM) {
    ggcm_mhd_calc_divb_bgrid_fc_ggcm(mhd, f, d);
  } else {
    assert(0);
  }
  
#if 0
  // FIXME, do we want to keep this?
  // If so, needs mrc_fld_mul() (pointwise multiplication)
  if (mhd->ymask) {
    struct mrc_fld *ymask = mrc_fld_get_as(mhd->ymask, FLD_TYPE);
    mrc_fld_mul(divb, ymask);
    mrc_fld_put_as(ymask, mhd->ymask);
  }
#endif

  double max_divb = mrc_fld_norm(d);
  mpi_printf(ggcm_mhd_comm(mhd), "max divb = %g\n", max_divb);

  mrc_fld_put_as(f, fld);
  mrc_fld_put_as(d, divb);

  if (0) {
    pde_free();
  }
}

