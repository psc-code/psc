
#include "ggcm_mhd_defs.h"
#include "ggcm_mhd_private.h"
#include "ggcm_mhd_crds.h"

#include <mrc_fld.h>
#include <mrc_domain.h>

#include <stdio.h>
#include <math.h>
#include <assert.h>

#include <mrc_fld_as_double.h>

// ----------------------------------------------------------------------
// ggcm_mhd_calc_divb

void
ggcm_mhd_calc_divb(struct ggcm_mhd *mhd, struct mrc_fld *fld, struct mrc_fld *divb)
{
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
}

