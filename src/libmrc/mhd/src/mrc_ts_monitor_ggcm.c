
#include <mrc_ts_monitor_private.h>

#include "ggcm_mhd_private.h"
#include "ggcm_mhd_commu.h"
#include "ggcm_mhd_diag.h"

#include <mrc_ts_private.h> // FIXME?
#include <mrc_io.h>

// ======================================================================
// mrc_ts_monitor subclass "ggcm"

struct mrc_ts_monitor_ggcm {
  int nr;
};

#define mrc_ts_monitor_ggcm(mon) mrc_to_subobj(mon, struct mrc_ts_monitor_ggcm)

static void
mrc_ts_monitor_ggcm_run(struct mrc_ts_monitor *mon, struct mrc_ts *ts)
{
  struct mrc_ts_monitor_ggcm *out = mrc_ts_monitor_ggcm(mon);

  mpi_printf(mrc_ts_monitor_comm(mon), "Writing output %d (time = %g)\n",
	     out->nr, ts->time);
  struct ggcm_mhd *mhd = (struct ggcm_mhd *) ts->ctx_obj;

  // FIXME? hack misusing istep as output counter
  mhd->istep = out->nr++;
  mhd->time = ts->time;
  ggcm_mhd_diag_run(mhd->diag);
}

struct mrc_ts_monitor_ops mrc_ts_monitor_ggcm_ops = {
  .name             = "ggcm",
  .size             = sizeof(struct mrc_ts_monitor_ggcm),
  .run              = mrc_ts_monitor_ggcm_run,
};

