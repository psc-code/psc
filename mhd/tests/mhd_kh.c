
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
#include <mrc_crds.h> 

#include <math.h>
#include <string.h>
#include <stdio.h>
#include <stdlib.h> 
#include <assert.h>

// ======================================================================
// ggcm_mhd_ic subclass "kh"

struct ggcm_mhd_ic_kh {
  float pert; // initial pertubation amplitude
  float r0; // initial density 0 
  float r1; // initial density 1  
  float v0; // velocity 0 
  float v1; // velocity 1
  float B0; // initial B
  float p0; // initial pressure
  float lambda; // wave number 
};

// ----------------------------------------------------------------------
// ggcm_mhd_ic_kh_run

static void
ggcm_mhd_ic_kh_run(struct ggcm_mhd_ic *ic)
{
  struct ggcm_mhd *mhd = ic->mhd;
  struct mrc_crds *crds = mrc_domain_get_crds(mhd->domain);  
  struct mrc_fld *fld = mrc_fld_get_as(mhd->fld, "mhd_pr_float");
  struct ggcm_mhd_ic_kh *sub = mrc_to_subobj(ic, struct ggcm_mhd_ic_kh);
  // FIXME, the "1" no of ghosts is ugly here, and caused by the use of
  // the B1* macros which shift the index (due to staggering)...
  // FIXME, need to set all components, can't rely on things being initialized to
  // zero because of the -> primitive conversion which divides by RR :(
  float xl[3], xh[3],  xmid[3], L[3], r[3];
  mrc_crds_get_param_float3(crds, "l", xl);
  mrc_crds_get_param_float3(crds, "h", xh);
  for(int i=0; i<3; i++){
    L[i] = xh[i] - xl[i];
    xmid[i] = 0.5 * (xh[i] + xl[i]);
  }
  //float wave1; 

  //float tmp = sub-> lambda; 
  mrc_fld_foreach(fld, ix,iy,iz, 1, 1) {
    r[0] = MRC_CRD(crds, 0, ix); 
    r[1] = MRC_CRD(crds, 1, iy);

    float wave1 = sin(2.* sub->lambda * M_PI * (r[0] - xmid[0]) / L[0]); 

    if(fabs(r[1]) < 0.25*L[1]){
      RR1(fld, ix,iy,iz) = sub->r0;
      PP1(fld, ix,iy,iz) = sub->p0;
      V1X(fld, ix,iy,iz) = sub->v0;
      V1Y(fld, ix,iy,iz) = sub->pert*wave1; 
    }else{
      RR1(fld, ix,iy,iz) = sub->r1;
      PP1(fld, ix,iy,iz) = sub->p0;
      V1X(fld, ix,iy,iz) = sub->v1;
      V1Y(fld, ix,iy,iz) = sub->pert*wave1; 
    }   
    B1X(fld, ix,iy,iz) = sub->B0; 
    B1Y(fld, ix,iy,iz) = 0.0; 
    V1Z(fld, ix,iy,iz) = 0.0;
  } mrc_fld_foreach_end;

  mrc_fld_put_as(fld, mhd->fld);
}



// ----------------------------------------------------------------------
// ggcm_mhd_ic_kh_descr

#define VAR(x) (void *)offsetof(struct ggcm_mhd_ic_kh, x)
static struct param ggcm_mhd_ic_kh_descr[] = {
  {"pert", VAR(pert), PARAM_FLOAT(1e-2)},  
  {"r0", VAR(r0), PARAM_FLOAT(2.0)},
  {"r1", VAR(r1), PARAM_FLOAT(1.0)},
  {"v0", VAR(v0), PARAM_FLOAT(0.5)},
  {"v1", VAR(v1), PARAM_FLOAT(-0.5)},
  {"B0", VAR(B0), PARAM_FLOAT(0.0)}, 
  {"p0", VAR(p0), PARAM_FLOAT(2.5)},
  {"lambda", VAR(lambda), PARAM_FLOAT(1.0)},
  {},
};
#undef VAR

// ----------------------------------------------------------------------
// ggcm_mhd_ic_kh_ops

struct ggcm_mhd_ic_ops ggcm_mhd_ic_kh_ops = {
  .name        = "kh",
  .size        = sizeof(struct ggcm_mhd_ic_kh),
  .param_descr = ggcm_mhd_ic_kh_descr,
  .run         = ggcm_mhd_ic_kh_run,
};



// ======================================================================
// ggcm_mhd class "kh"

// ----------------------------------------------------------------------
// ggcm_mhd_kh_create

static void
ggcm_mhd_kh_create(struct ggcm_mhd *mhd)
{
  ggcm_mhd_default_box(mhd);

  /* set defaults for coord arrays */
  struct mrc_crds *crds = mrc_domain_get_crds(mhd->domain);
  mrc_crds_set_type(crds, "uniform");
  mrc_crds_set_param_int(crds, "sw", SW_2);   // 'stencil width' 
  mrc_crds_set_param_float3(crds, "l", (float[3]) {  0.0, 0.0, 0.0 });
  mrc_crds_set_param_float3(crds, "h", (float[3]) {  1.0, 1.0, 0.1 });
}

static struct ggcm_mhd_ops ggcm_mhd_kh_ops = {
  .name             = "kh",
  .create           = ggcm_mhd_kh_create,
};

// ======================================================================

extern struct ggcm_mhd_diag_ops ggcm_mhd_diag_c_ops;

int
main(int argc, char **argv)
{
  MPI_Init(&argc, &argv);
  libmrc_params_init(argc, argv);
  ggcm_mhd_register();

  mrc_class_register_subclass(&mrc_class_ggcm_mhd, &ggcm_mhd_kh_ops);  
  mrc_class_register_subclass(&mrc_class_ggcm_mhd_diag, &ggcm_mhd_diag_c_ops);

  mrc_class_register_subclass(&mrc_class_ggcm_mhd_ic, &ggcm_mhd_ic_kh_ops);  
 
  struct ggcm_mhd *mhd = ggcm_mhd_create(MPI_COMM_WORLD);
  ggcm_mhd_set_type(mhd, "kh");
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

