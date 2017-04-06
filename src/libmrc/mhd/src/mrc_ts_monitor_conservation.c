
#include <mrc_ts_monitor_private.h>

#include "ggcm_mhd_private.h"
#include "ggcm_mhd_diag.h"
#include "ggcm_mhd_defs.h"

#include <mrc_ts_private.h> // FIXME?
#include <mrc_io.h>
#include <mrc_fld_as_double.h>

// ----------------------------------------------------------------------
// ggcm_mhd_calc_integrals

static void
ggcm_mhd_calc_integrals(struct ggcm_mhd *mhd, double *vals)
{
  int gdims[3]; mrc_domain_get_global_dims(mhd->domain, gdims);
  struct mrc_crds *crds = mrc_domain_get_crds(mhd->domain);

  struct mrc_fld *f = mrc_fld_get_as(mhd->fld, FLD_TYPE);

  for (int m = 0; m < 5; m++) {
    vals[m] = 0.;
  }
  for (int p = 0; p < mrc_fld_nr_patches(f); p++) {
    double dx[3]; mrc_crds_get_dx(crds, p, dx);
    double dV = 1.;
    for (int d = 0; d < 3; d++) {
      if (gdims[d] > 1) {
	dV *= dx[d];
      }
    }
    mrc_fld_foreach(f, ix,iy,iz, 0, 0) {
      for (int m = 0; m < 5; m++) {
	vals[m] += dV * M3(f, m, ix,iy,iz, p);
      }
    } mrc_fld_foreach_end;
  }

  mrc_fld_put_as(f, mhd->fld);
}

// ======================================================================
// mrc_ts_monitor subclass "conservation"

struct mrc_ts_monitor_conservation {
  FILE *file;
};

#define mrc_ts_monitor_conservation(mon) mrc_to_subobj(mon, struct mrc_ts_monitor_conservation)

// ----------------------------------------------------------------------
// mrc_ts_monitor_conservation_setup

static void
mrc_ts_monitor_conservation_setup(struct mrc_ts_monitor *mon)
{
  struct mrc_ts_monitor_conservation *sub = mrc_ts_monitor_conservation(mon);

  int size; MPI_Comm_size(mrc_ts_monitor_comm(mon), &size);
  assert(size == 1);
  
  sub->file = fopen("diag.asc", "w");
}

// ----------------------------------------------------------------------
// mrc_ts_monitor_conservation_destroy

static void
mrc_ts_monitor_conservation_destroy(struct mrc_ts_monitor *mon)
{
  struct mrc_ts_monitor_conservation *sub = mrc_ts_monitor_conservation(mon);

  if (sub->file) {
    fclose(sub->file);
  }
}

// ----------------------------------------------------------------------
// mrc_ts_monitor_conservation_run

static void
mrc_ts_monitor_conservation_run(struct mrc_ts_monitor *mon, struct mrc_ts *ts)
{
  struct mrc_ts_monitor_conservation *sub = mrc_ts_monitor_conservation(mon);
  struct ggcm_mhd *mhd = (struct ggcm_mhd *) ts->ctx_obj;

  mpi_printf(mrc_ts_monitor_comm(mon), "Monitoring conservation (time = %g)\n",
	     ts->time);

  mhd->time_code = ts->time;
  ggcm_mhd_fill_ghosts(mhd, mhd->fld, mhd->time_code);

  double vals[5];
  ggcm_mhd_calc_integrals(mhd, vals);
  
  assert(sub->file);
  fprintf(sub->file, "%g", ts->time);
  for (int m = 0; m < 5; m++) {
    fprintf(sub->file, " %g", vals[m]);
  }
  fprintf(sub->file, "\n");
  fflush(sub->file);
}

struct mrc_ts_monitor_ops mrc_ts_monitor_conservation_ops = {
  .name             = "conservation",
  .size             = sizeof(struct mrc_ts_monitor_conservation),
  .setup            = mrc_ts_monitor_conservation_setup,
  .destroy          = mrc_ts_monitor_conservation_destroy,
  .run              = mrc_ts_monitor_conservation_run,
};

