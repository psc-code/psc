
//#define BOUNDS_CHECK

#include "ggcm_mhd_defs.h"
#include "ggcm_mhd_private.h"
#include "ggcm_mhd_step.h"
#include "ggcm_mhd_crds_private.h"
#include "ggcm_mhd_crds_gen.h"
#include "ggcm_mhd_bnd.h"
#include "ggcm_mhd_diag.h"
#include "ggcm_mhd_ic_private.h"

#include <mrc_ts.h>
#include <mrc_ts_monitor.h>
#include <mrc_fld.h>
#include <mrc_domain.h>
#include <mrc_params.h>
#include <mrc_ddc.h>
#include <mrctest.h>
#include <mrc_io.h> 

#include <math.h>
#include <string.h>
#include <stdio.h>
#include <stdlib.h> 
#include <assert.h>

struct ggcm_mhd_ic_whistler {
  float mpermi;
  float Boz;
  float pert; 
  float eps; 
  float n0; 
  float lambda;   
};
// ----------------------------------------------------------------------
// ggcm_mhd_ic_whistler_run

static void
ggcm_mhd_ic_whistler_run(struct ggcm_mhd_ic *ic)
{
  struct ggcm_mhd_ic_whistler *sub = mrc_to_subobj(ic, struct ggcm_mhd_ic_whistler);
  struct ggcm_mhd *gmhd = ic->mhd;  
  struct mrc_fld *f3 = mrc_fld_get_as(gmhd->fld, "float");
  struct mrc_crds *crds = mrc_domain_get_crds(gmhd->domain);  
  float xl[3], xh[3], L[3], r[3];
  mrc_crds_get_param_float3(crds, "l", xl);
  mrc_crds_get_param_float3(crds, "h", xh);
  for(int i=0; i<3; i++){
    L[i] = xh[i] - xl[i];
  }
  float gamma = gmhd->par.gamm;
  mrc_fld_foreach(f3, ix, iy, iz, 1, 1) {
    r[0] = MRC_CRD(crds, 0, ix);
    r[1] = MRC_CRD(crds, 1, iy);
    r[2] = MRC_CRD(crds, 2, iz);

    // Xin etal (2005)      
    float kk= (sub->lambda * 2.*M_PI) / L[2] ;
    const int *dims = mrc_fld_dims(f3);
    int nz = dims[2];
    float vp= 2.*M_PI*sub->Boz*nz/((sub->mpermi)*(sub->n0)*L[2]); 
    MRC_F3(f3, _B1X, ix+1,iy,iz) = (sub->pert) * vp * sin( kk*r[2] ) ;       
    MRC_F3(f3, _B1Y, ix,iy+1,iz) = -(sub->pert) * vp * cos( kk*r[2] ) ;   
    MRC_F3(f3, _B1Z, ix,iy,iz+1) = sub->Boz ; 
    MRC_F3(f3, _RV1X, ix,iy,iz) = (sub->pert) * sub->n0 * sin( kk*r[2] ) ;      
    MRC_F3(f3, _RV1Y, ix,iy,iz) = -(sub->pert) * sub->n0 * cos( kk*r[2] ) ;
    MRC_F3(f3, _RR1, ix, iy, iz) = sub->n0 ;       

    MRC_F3(f3, _UU1 , ix, iy, iz) =  MRC_F3(f3, _RR1, ix, iy, iz) / (gamma -1.f) + 	
      .5f * (sqr(MRC_F3(f3, _RV1X, ix, iy, iz)) +
	     sqr(MRC_F3(f3, _RV1Y, ix, iy, iz)) +
	     sqr(MRC_F3(f3, _RV1Z, ix, iy, iz))) / MRC_F3(f3, _RR1, ix, iy, iz) +
      .5f * (sqr(.5*(B1X(f3, ix,iy,iz) + B1X(f3, ix+1,iy,iz))) +
	     sqr(.5*(B1Y(f3, ix,iy,iz) + B1Y(f3, ix,iy+1,iz))) +
	     sqr(.5*(B1Z(f3, ix,iy,iz) + B1Z(f3, ix,iy,iz+1))));
  } mrc_fld_foreach_end;
}

// ----------------------------------------------------------------------
// ggcm_mhd_ic_whistler_descr

#define VAR(x) (void *)offsetof(struct ggcm_mhd_ic_whistler, x)
static struct param ggcm_mhd_ic_whistler_descr[] = {
  {"mpermi", VAR(mpermi), PARAM_FLOAT(1.0)},
  {"pert", VAR(pert), PARAM_FLOAT(1e-5)},
  {"Boz", VAR(Boz), PARAM_FLOAT(1.0)},
  {"n0", VAR(n0), PARAM_FLOAT(25.)},
  {"lambda", VAR(lambda), PARAM_FLOAT(4)},  
  {},
};
#undef VAR



// ----------------------------------------------------------------------
// ggcm_mhd_ic_whistler_ops

struct ggcm_mhd_ic_ops ggcm_mhd_ic_whistler_ops = {
  .name        = "whistler",
  .size        = sizeof(struct ggcm_mhd_ic_whistler),
  .param_descr = ggcm_mhd_ic_whistler_descr,
  .run         = ggcm_mhd_ic_whistler_run,
};



// ======================================================================
// ggcm_mhd class "whistler"

// ----------------------------------------------------------------------
// ggcm_mhd_whistler_create

static void
ggcm_mhd_whistler_create(struct ggcm_mhd *mhd)
{
  ggcm_mhd_default_box(mhd);

  /* set defaults for coord arrays */
  struct mrc_crds *crds = mrc_domain_get_crds(mhd->domain);
  mrc_crds_set_type(crds, "uniform");
  mrc_crds_set_param_int(crds, "sw", SW_2);   // 'stencil width' 
  mrc_crds_set_param_float3(crds, "l", (float[3]) {  0.0, 0.0, 0.0 });
  mrc_crds_set_param_float3(crds, "h", (float[3]) {  2.0, 0.1, 0.1 });

  /* set defaults for the ddc, this does the communication */
  struct mrc_ddc *ddc = mrc_domain_get_ddc(mhd->domain);
  mrc_ddc_set_param_int(ddc, "max_n_fields", 8);
  mrc_ddc_set_param_int3(ddc, "ibn", (int[3]) { SW_2, SW_2, SW_2 });

  // generate MHD solver grid from mrc_crds
  ggcm_mhd_crds_gen_set_type(mhd->crds->crds_gen, "mrc");
}

static struct ggcm_mhd_ops ggcm_mhd_whistler_ops = {
  .name             = "whistler",
  .create           = ggcm_mhd_whistler_create,
};

// ======================================================================

extern struct ggcm_mhd_diag_ops ggcm_mhd_diag_c_ops;

int
main(int argc, char **argv)
{
  MPI_Init(&argc, &argv);
  libmrc_params_init(argc, argv);
  ggcm_mhd_register();

  mrc_class_register_subclass(&mrc_class_ggcm_mhd, &ggcm_mhd_whistler_ops);  
  mrc_class_register_subclass(&mrc_class_ggcm_mhd_diag, &ggcm_mhd_diag_c_ops);

  mrc_class_register_subclass(&mrc_class_ggcm_mhd_ic, &ggcm_mhd_ic_whistler_ops);  
 
  struct ggcm_mhd *mhd = ggcm_mhd_create(MPI_COMM_WORLD);
  ggcm_mhd_set_type(mhd, "whistler");
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

