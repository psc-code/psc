
//#define BOUNDS_CHECK

#include "ggcm_mhd_defs.h"
#include "ggcm_mhd_private.h"
#include "ggcm_mhd_step.h"
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

// ======================================================================
// ggcm_mhd_ic subclass "bw"

struct ggcm_mhd_ic_bw {
  float mpermi;
  const char* pdim;
};

// ----------------------------------------------------------------------
// ggcm_mhd_ic_bw_run

static void
ggcm_mhd_ic_bw_run(struct ggcm_mhd_ic *ic)
{
  struct ggcm_mhd_ic_bw *sub = mrc_to_subobj(ic, struct ggcm_mhd_ic_bw);
  struct ggcm_mhd *mhd = ic->mhd;
  //struct mrc_fld *f3 = mrc_fld_get_as(mhd->fld, "float");
  struct mrc_fld *fld = mrc_fld_get_as(mhd->fld, "mhd_pr_float");
  struct mrc_crds *crds = mrc_domain_get_crds(mhd->domain);  

  float xl[3], xh[3], L[3], r[3];
  mrc_crds_get_param_float3(crds, "l", xl);
  mrc_crds_get_param_float3(crds, "h", xh);
  for(int i = 0; i < 3; i++){
    L[i] = xh[i] - xl[i];
  }

  mrc_fld_foreach(fld, ix, iy, iz, 2, 2) {
    r[0] = MRC_CRD(crds, 0, ix);
    r[1] = MRC_CRD(crds, 1, iy);
    r[2] = MRC_CRD(crds, 2, iz);
    
    if(strcmp(sub->pdim, "x") == 0){

      MHERE;
      if(fabs(r[0]) < 0.5*L[0]){
	// Left                         
	RR1(fld, ix,iy,iz) = 1.0;
	PP1(fld, ix,iy,iz) = RR1(fld, ix,iy,iz);
	V1X(fld, ix,iy,iz) = 0.0;
	V1Y(fld, ix,iy,iz) = 0.0;
	V1Z(fld, ix,iy,iz) = 0.0;
	B1X(fld, ix,iy,iz) = 0.75; 
	B1Y(fld, ix,iy,iz) = 1.0;
	B1Z(fld, ix,iy,iz) = 0.0;
      } else {
	// Right
	// Left                         
	RR1(fld, ix,iy,iz) = 0.125;	 
	PP1(fld, ix,iy,iz) = RR1(fld, ix,iy,iz);
	V1X(fld, ix,iy,iz) = 0.0;
	V1Y(fld, ix,iy,iz) = 0.0;
	V1Z(fld, ix,iy,iz) = 0.0;
	B1X(fld, ix,iy,iz) = 0.75;
	B1Y(fld, ix,iy,iz) = -1.0;
	B1Z(fld, ix,iy,iz) = 0.0;


#if 0        

	/*
      if(fabs(r[0]) < 0.5*L[0]){
	// Left                         
	MRC_F3(f3, _RR1, ix, iy, iz) = 1.0;
	MRC_F3(f3, _B1X, ix, iy, iz) = 0.75;
	MRC_F3(f3, _B1Y, ix, iy, iz) = 1.0;
	MRC_F3(f3, _B1Z, ix, iy, iz) = 0.0;
	//B1X(f3, ix,iy,iz) = 0.75;
	//B1Y(f3, ix,iy,iz) = 1.0;
	//B1Z(f3, ix,iy,iz) = 0.0; 
	
	MRC_F3(f3, _UU1 , ix, iy, iz) = 1.0 / (gamma - 1.f) +
	  .5f * (sqr(MRC_F3(f3, _RV1X, ix, iy, iz)) +
		 sqr(MRC_F3(f3, _RV1Y, ix, iy, iz)) +
		 sqr(MRC_F3(f3, _RV1Z, ix, iy, iz))) / MRC_F3(f3, _RR1, ix, iy, iz) +
	  .5f * (sqr(.5*(B1X(f3, ix,iy,iz) + B1X(f3, ix+1,iy,iz))) +
		 sqr(.5*(B1Y(f3, ix,iy,iz) + B1Y(f3, ix,iy+1,iz))) +
		 sqr(.5*(B1Z(f3, ix,iy,iz) + B1Z(f3, ix,iy,iz+1))));
	
      } else {
	// Right
	MRC_F3(f3, _RR1, ix, iy, iz) = 0.125;
	MRC_F3(f3, _B1X, ix, iy, iz) = 0.75;
	MRC_F3(f3, _B1Y, ix, iy, iz) = -1.0;
	MRC_F3(f3, _B1Z, ix, iy, iz) = 0.0;

	MRC_F3(f3, _UU1 , ix, iy, iz) = 0.1 / (gamma - 1.f) +
	  .5f * (sqr(MRC_F3(f3, _RV1X, ix, iy, iz)) +
		 sqr(MRC_F3(f3, _RV1Y, ix, iy, iz)) +
		 sqr(MRC_F3(f3, _RV1Z, ix, iy, iz))) / MRC_F3(f3, _RR1, ix, iy, iz) +
	  .5f * (sqr(.5*(B1X(f3, ix,iy,iz) + B1X(f3, ix+1,iy,iz))) +
		 sqr(.5*(B1Y(f3, ix,iy,iz) + B1Y(f3, ix,iy+1,iz))) +
		 sqr(.5*(B1Z(f3, ix,iy,iz) + B1Z(f3, ix,iy,iz+1))));
	*/
#endif


      }
  } else if(strcmp(sub->pdim, "y") == 1){

      /*
    if(fabs(r[1]) < 0.5*L[1]){
      // Left 
      MRC_F3(f3, _RR1, ix, iy, iz) = 1.0;
      MRC_F3(f3, _B1Y , ix, iy, iz) = 0.75;
      MRC_F3(f3, _B1X , ix, iy, iz) = 1.0;
      MRC_F3(f3, _B1Z , ix, iy, iz) = 0.0;
      MRC_F3(f3, _UU1 , ix, iy, iz) = 1.0 / (gamma - 1.f) +
	.5f * (sqr(MRC_F3(f3, _RV1X, ix, iy, iz)) +
	       sqr(MRC_F3(f3, _RV1Y, ix, iy, iz)) +
	       sqr(MRC_F3(f3, _RV1Z, ix, iy, iz))) / MRC_F3(f3, _RR1, ix, iy, iz) +
	.5f * (sqr(MRC_F3(f3, _B1X, ix, iy, iz)) +
	       sqr(MRC_F3(f3, _B1Y, ix, iy, iz)) +
	       sqr(MRC_F3(f3, _B1Z, ix, iy, iz)));      
    }else{
      // Right
      MRC_F3(f3, _RR1, ix, iy, iz) = 0.125;
      MRC_F3(f3, _B1Y , ix, iy, iz) = 0.75;
      MRC_F3(f3, _B1X , ix, iy, iz) = -1.0;
      MRC_F3(f3, _B1Z , ix, iy, iz) = 0.0;
      MRC_F3(f3, _UU1 , ix, iy, iz) = 0.1 / (gamma - 1.f) +
	.5f * (sqr(MRC_F3(f3, _RV1X, ix, iy, iz)) +
	       sqr(MRC_F3(f3, _RV1Y, ix, iy, iz)) +
	       sqr(MRC_F3(f3, _RV1Z, ix, iy, iz))) / MRC_F3(f3, _RR1, ix, iy, iz) +
	.5f * (sqr(MRC_F3(f3, _B1X, ix, iy, iz)) +
	       sqr(MRC_F3(f3, _B1Y, ix, iy, iz)) +
	       sqr(MRC_F3(f3, _B1Z, ix, iy, iz)));      
    }
    MRC_F3(f3, _RV1Z , ix, iy, iz) = 0.0;
    MRC_F3(f3, _RV1X , ix, iy, iz) = 0.0;
    MRC_F3(f3, _RV1Y , ix, iy, iz) = 0.0;
    
  } else if(strcmp(sub->pdim, "z") == 1){
    if(fabs(r[2]) < 0.5*L[2]){
      // Left 
      MRC_F3(f3, _RR1, ix, iy, iz) = 1.0;
      MRC_F3(f3, _RV1X , ix, iy, iz) = 0.0;
      MRC_F3(f3, _RV1Y , ix, iy, iz) = 0.0;	
      MRC_F3(f3, _B1Y , ix, iy, iz) = 0.0;
      MRC_F3(f3, _B1Z , ix, iy, iz) = 0.75;
      MRC_F3(f3, _B1X , ix, iy, iz) = 1.0;
      MRC_F3(f3, _UU1 , ix, iy, iz) = 1.0 / (gamma - 1.f) +
	.5f * (sqr(MRC_F3(f3, _RV1X, ix, iy, iz)) +
	       sqr(MRC_F3(f3, _RV1Y, ix, iy, iz)) +
	       sqr(MRC_F3(f3, _RV1Z, ix, iy, iz))) / MRC_F3(f3, _RR1, ix, iy, iz) +
	.5f * (sqr(MRC_F3(f3, _B1X, ix, iy, iz)) +
	       sqr(MRC_F3(f3, _B1Y, ix, iy, iz)) +
	       sqr(MRC_F3(f3, _B1Z, ix, iy, iz)));      
    }else{
      // Right
      MRC_F3(f3, _RR1, ix, iy, iz) = 0.125;
      MRC_F3(f3, _B1Y , ix, iy, iz) = 0.0;
      MRC_F3(f3, _B1X , ix, iy, iz) = -1.0;
      MRC_F3(f3, _B1Z , ix, iy, iz) = 0.75;
      MRC_F3(f3, _UU1 , ix, iy, iz) = 0.1 / (gamma - 1.f) +
	.5f * (sqr(MRC_F3(f3, _RV1X, ix, iy, iz)) +
	       sqr(MRC_F3(f3, _RV1Y, ix, iy, iz)) +
	       sqr(MRC_F3(f3, _RV1Z, ix, iy, iz))) / MRC_F3(f3, _RR1, ix, iy, iz) +
	.5f * (sqr(MRC_F3(f3, _B1X, ix, iy, iz)) +
	       sqr(MRC_F3(f3, _B1Y, ix, iy, iz)) +
	       sqr(MRC_F3(f3, _B1Z, ix, iy, iz)));      
    }
      */

      /*
    MRC_F3(f3, _RV1X , ix, iy, iz) = 0.0;
    MRC_F3(f3, _RV1Y , ix, iy, iz) = 0.0;	
    MRC_F3(f3, _RV1Z , ix, iy, iz) = 0.0;       
      */

  } else {           
    assert(0); /* unknown initial condition */
  }
  } mrc_fld_foreach_end;
  
  mrc_fld_put_as(fld, mhd->fld);
}

// ----------------------------------------------------------------------
// ggcm_mhd_ic_bw_descr

#define VAR(x) (void *)offsetof(struct ggcm_mhd_ic_bw, x)
static struct param ggcm_mhd_ic_bw_descr[] = {
  {"mpermi", VAR(mpermi), PARAM_FLOAT(1.)},
  {"pdim", VAR(pdim), PARAM_STRING("x")},  
  {},
};
#undef VAR

// ----------------------------------------------------------------------
// ggcm_mhd_ic_bw_ops

struct ggcm_mhd_ic_ops ggcm_mhd_ic_bw_ops = {
  .name        = "bw",
  .size        = sizeof(struct ggcm_mhd_ic_bw),
  .param_descr = ggcm_mhd_ic_bw_descr,
  .run         = ggcm_mhd_ic_bw_run,
};


// ======================================================================
// ggcm_mhd class "bw"

// ----------------------------------------------------------------------
// ggcm_mhd_bw_create

static void
ggcm_mhd_bw_create(struct ggcm_mhd *mhd)
{
  ggcm_mhd_default_box(mhd);

  ggcm_mhd_bnd_set_type(mhd->bnd, "open_x");
  mrc_domain_set_param_int(mhd->domain, "bcx", BC_NONE);

  /* set defaults for coord arrays */
  struct mrc_crds *crds = mrc_domain_get_crds(mhd->domain);
  mrc_crds_set_type(crds, "uniform");
  mrc_crds_set_param_int(crds, "sw", SW_2);   // 'stencil width' 
  mrc_crds_set_param_float3(crds, "l", (float[3]) {  0.0, 0.0, 0.0 });
  mrc_crds_set_param_float3(crds, "h", (float[3]) {  1.0, 0.05, 0.05 });

  /* set defaults for the ddc, this does the communication */
  struct mrc_ddc *ddc = mrc_domain_get_ddc(mhd->domain);
  mrc_ddc_set_param_int(ddc, "max_n_fields", 8);
  mrc_ddc_set_param_int3(ddc, "ibn", (int[3]) { SW_2, SW_2, SW_2 });
}

static struct ggcm_mhd_ops ggcm_mhd_bw_ops = {
  .name             = "bw",
  .create           = ggcm_mhd_bw_create,
};

// ======================================================================

extern struct ggcm_mhd_diag_ops ggcm_mhd_diag_c_ops;

int
main(int argc, char **argv)
{
  MPI_Init(&argc, &argv);
  libmrc_params_init(argc, argv);
  ggcm_mhd_register();

  mrc_class_register_subclass(&mrc_class_ggcm_mhd, &ggcm_mhd_bw_ops);  
  mrc_class_register_subclass(&mrc_class_ggcm_mhd_diag, &ggcm_mhd_diag_c_ops);
  mrc_class_register_subclass(&mrc_class_ggcm_mhd_ic, &ggcm_mhd_ic_bw_ops);  
 
  struct ggcm_mhd *mhd = ggcm_mhd_create(MPI_COMM_WORLD);
  ggcm_mhd_set_type(mhd, "bw");
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

