
//#define BOUNDS_CHECK

#include "ggcm_mhd_defs.h"
#include "ggcm_mhd_private.h"
#include "ggcm_mhd_step.h"
#include "ggcm_mhd_crds.h"
#include "ggcm_mhd_diag.h"
#include "ggcm_mhd_ic_private.h"

#include <mrc_ts.h>
#include <mrc_ts_monitor.h>
#include <mrc_fld.h>
#include <mrc_domain.h>
#include <mrc_ddc.h>

#include <math.h>
#include <string.h>
#include <stdio.h>
#include <stdlib.h> 
#include <assert.h>

// ======================================================================
// ggcm_mhd_ic subclass "ot"

// ----------------------------------------------------------------------
// ggcm_mhd_ic_ot_run

static void
ggcm_mhd_ic_ot_run(struct ggcm_mhd_ic *ic)
{
  //  struct ggcm_mhd_ic_ot *sub = mrc_to_subobj(ic, struct ggcm_mhd_ic_ot);

  struct ggcm_mhd *mhd = ic->mhd;
  struct mrc_crds *crds = mrc_domain_get_crds(mhd->domain);  
  struct mrc_fld *fld = mrc_fld_get_as(mhd->fld, "mhd_pr_float");

  // FIXME, the "1" no of ghosts is ugly here, and caused by the use of
  // the B1* macros which shift the index (due to staggering)...
  // FIXME, need to set all components, can't rely on things being initialized to
  // zero because of the -> primitive conversion which divides by RR :(
  mrc_fld_foreach(fld, ix,iy,iz, 1, 1) {
    float r[3];
    r[0] = MRC_CRD(crds, 0, ix);
    r[1] = MRC_CRD(crds, 1, iy);
    RR1(fld, ix,iy,iz) = 25. / (36.*M_PI);
    PP1(fld, ix,iy,iz) = RR1(fld, ix,iy,iz);
    V1X(fld, ix,iy,iz) = - sin(2. * M_PI * r[1]);
    V1Y(fld, ix,iy,iz) =   sin(2. * M_PI * r[0]);
    V1Z(fld, ix,iy,iz) = 0.;
    B1X(fld, ix,iy,iz) = - sqrt(1./(4.*M_PI)) * sin(2. * M_PI * r[1]); 
    B1Y(fld, ix,iy,iz) =   sqrt(1./(4.*M_PI)) * sin(4. * M_PI * r[0]);
    B1Z(fld, ix,iy,iz) = 0.;
  } mrc_fld_foreach_end;

  mrc_fld_put_as(fld, mhd->fld);
}

// ----------------------------------------------------------------------
// ggcm_mhd_ic_ot_ops

struct ggcm_mhd_ic_ops ggcm_mhd_ic_ot_ops = {
  .name        = "ot",
  .run         = ggcm_mhd_ic_ot_run,
};



// ======================================================================
// ggcm_mhd class "ot"

// ----------------------------------------------------------------------
// ggcm_mhd_ot_create

static void
ggcm_mhd_ot_create(struct ggcm_mhd *mhd)
{
  ggcm_mhd_default_box(mhd);

  /* set defaults for coord arrays */
  struct mrc_crds *crds = mrc_domain_get_crds(mhd->domain);
  mrc_crds_set_type(crds, "uniform");
  mrc_crds_set_param_int(crds, "sw", SW_2);   // 'stencil width' 
  mrc_crds_set_param_float3(crds, "l", (float[3]) {  0.0, 0.0, 0.0 });
  mrc_crds_set_param_float3(crds, "h", (float[3]) {  1.0, 1.0, 0.1 });

  /* set defaults for the ddc, this does the communication */
  struct mrc_ddc *ddc = mrc_domain_get_ddc(mhd->domain);
  mrc_ddc_set_param_int(ddc, "max_n_fields", 8);
  mrc_ddc_set_param_int3(ddc, "ibn", (int[3]) { SW_2, SW_2, SW_2 });
}

static struct ggcm_mhd_ops ggcm_mhd_ot_ops = {
  .name             = "ot",
  .create           = ggcm_mhd_ot_create,
};

// ======================================================================

extern struct ggcm_mhd_diag_ops ggcm_mhd_diag_c_ops;

int
main(int argc, char **argv)
{
  MPI_Init(&argc, &argv);
  libmrc_params_init(argc, argv);
  ggcm_mhd_register();

  mrc_class_register_subclass(&mrc_class_ggcm_mhd, &ggcm_mhd_ot_ops);  
  mrc_class_register_subclass(&mrc_class_ggcm_mhd_diag, &ggcm_mhd_diag_c_ops);

  mrc_class_register_subclass(&mrc_class_ggcm_mhd_ic, &ggcm_mhd_ic_ot_ops);  
 
  struct ggcm_mhd *mhd = ggcm_mhd_create(MPI_COMM_WORLD);
  ggcm_mhd_set_type(mhd, "ot");
  mrc_fld_set_type(mhd->fld, "mhd_fc_float");
  ggcm_mhd_step_set_type(mhd->step, "cweno");
  ggcm_mhd_set_from_options(mhd);
  ggcm_mhd_setup(mhd);
  ggcm_mhd_view(mhd);

  // set up initial condition
  ggcm_mhd_ic_run(mhd->ic);

  // run time integration
  struct mrc_ts *ts = mrc_ts_create(mrc_domain_comm(mhd->domain));
  mrc_ts_set_type(ts, "rk2");
  mrc_ts_set_context(ts, ggcm_mhd_to_mrc_obj(mhd));

  struct mrc_ts_monitor *mon_output =
    mrc_ts_monitor_create(mrc_ts_comm(ts));
  mrc_ts_monitor_set_type(mon_output, "ggcm");
  mrc_ts_monitor_set_name(mon_output, "mrc_ts_output");
  mrc_ts_add_monitor(ts, mon_output);

  mrc_ts_set_dt(ts, 1e-6);
  mrc_ts_set_solution(ts, mrc_fld_to_mrc_obj(mhd->fld));
  mrc_ts_set_rhs_function(ts, ts_ggcm_mhd_step_calc_rhs, mhd);
  mrc_ts_set_from_options(ts);
  mrc_ts_view(ts);
  mrc_ts_setup(ts);
  mrc_ts_solve(ts);
  mrc_ts_view(ts);
  mrc_ts_destroy(ts);  
  ggcm_mhd_destroy(mhd);

  MPI_Finalize();
  return 0;
}

