
#include "ggcm_mhd_defs.h"
#include "ggcm_mhd_private.h"
#include "ggcm_mhd_bnd.h"
#include "ggcm_mhd_step.h"

#include <mrc_ts_private.h>
#include <mrc_profile.h>

// ======================================================================
// mrc_ts subclass "ggcm"

// ----------------------------------------------------------------------
// mrc_ts_ggcm_step

static void
mrc_ts_ggcm_step(struct mrc_ts *ts)
{
  static int PR_push;
  if (!PR_push) {
    PR_push = prof_register("mrc_ts_ggcm_step", 1., 0, 0);
  }

  prof_start(PR_push);

  struct ggcm_mhd *mhd = (struct ggcm_mhd *) ts->ctx_obj;
  struct mrc_fld *fld = (struct mrc_fld *) ts->x;
  assert(fld == mhd->fld);
  mhd->time = ts->time;
  mhd->dt = ts->dt;

  ggcm_mhd_fill_ghosts(mhd, fld, _RR1, mhd->time);
  ggcm_mhd_step_pred(mhd->step);

  ggcm_mhd_fill_ghosts(mhd, fld, _RR2, mhd->time + mhd->bndt);
  ggcm_mhd_step_corr(mhd->step);

  prof_stop(PR_push);
}

// ----------------------------------------------------------------------
// mrc_ts_ggcm_ops

struct mrc_ts_ops mrc_ts_ggcm_ops = {
  .name             = "ggcm",
  .step             = mrc_ts_ggcm_step,
};

