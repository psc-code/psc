
//#define BOUNDS_CHECK

#include "ggcm_mhd_defs.h"
#include "ggcm_mhd_private.h"
#include "ggcm_mhd_step.h"
#include "ggcm_mhd_crds.h"
#include "ggcm_mhd_bnd.h"
#include "ggcm_mhd_diag.h"
#include "ggcm_mhd_ic_private.h"

#include <mrc_ts.h>
#include <mrc_ts_monitor.h>
#include <mrc_fld.h>
#include <mrc_domain.h>

#include <math.h>
#include <string.h>
#include <stdio.h>
#include <stdlib.h> 
#include <assert.h>

// ======================================================================
// ggcm_mhd_ic subclass "fadeev"

struct ggcm_mhd_ic_fadeev {
  float Bo; 
  float Boz;
  float pert; 
  float eps; 
  float lambda;   
  float dens0; 
};
// ----------------------------------------------------------------------
// ggcm_mhd_ic_fadeev_run

static void
ggcm_mhd_ic_fadeev_run(struct ggcm_mhd_ic *ic)
{
  struct ggcm_mhd_ic_fadeev *sub = mrc_to_subobj(ic, struct ggcm_mhd_ic_fadeev);
  struct ggcm_mhd *mhd = ic->mhd;
  struct mrc_crds *crds = mrc_domain_get_crds(mhd->domain);  

  struct mrc_fld *fld_psi = mrc_domain_fld_create(mhd->domain, SW_2, "psi");
  mrc_fld_setup(fld_psi);

  float xl[3], xh[3], L[3], r[3];
  mrc_crds_get_param_float3(crds, "l", xl);
  mrc_crds_get_param_float3(crds, "h", xh);
  for(int i = 0; i < 3; i++){
    L[i] = xh[i] - xl[i];
  }
  
  float Bo = sub->Bo;
  float pert = sub->pert;
  float Boz = sub->Boz;
  float eps = sub->eps;
  float lam = (sub->lambda)*L[0] ;  // defines island size   
  float kk = (2.*M_PI) / lam ;      

  struct mrc_fld *psi = mrc_fld_get_as(fld_psi, "float");
  struct mrc_fld *fld = mrc_fld_get_as(mhd->fld, "mhd_pr_float");

  mrc_fld_foreach(psi, ix, iy, iz, 1, 2) {
    r[0] = .5*(MRC_CRDX(crds, ix) + MRC_CRDX(crds, ix-1));
    r[1] = .5*(MRC_CRDY(crds, iy) + MRC_CRDY(crds, iy-1));
    
    MRC_F3(psi, 0, ix,iy,iz) = -(Bo / kk)*( log(cosh(kk*r[1]) + eps*cos(kk*r[0])));
  } mrc_fld_foreach_end;

  float *bd2x = ggcm_mhd_crds_get_crd(mhd->crds, 0, BD2);
  float *bd2y = ggcm_mhd_crds_get_crd(mhd->crds, 1, BD2);

  mrc_fld_foreach(fld, ix, iy, iz, 1, 1) {
    // FIXME! the staggering for B is okay, but fld_psi and other stuff below needs to be
    // fixed / checked for cell-centered
    r[0] = MRC_CRD(crds, 0, ix);
    r[1] = MRC_CRD(crds, 1, iy);

    RR1(fld, ix,iy,iz)  = 0.5*sqr(Bo) * (1.0-sqr(eps)) * 
      exp(2.0*kk* MRC_F3(fld_psi, 0, ix,iy,iz)/(Bo)) + 0.5*sqr(Boz) + sub->dens0;
    PP1(fld, ix,iy,iz) = RR1(fld, ix,iy,iz);
    V1X(fld, ix,iy,iz) = (pert) * (1.-kk*kk*r[0]*r[0]) *
      exp(-kk*kk*r[1]*r[1])*sin(kk*r[0]*0.5);	
    V1Y(fld, ix,iy,iz) = -(pert) * ( 0.5*kk*r[1] ) *
      exp(-kk*kk*r[1]*r[1])*cos(kk*r[0]*0.5);            
    V1Z(fld, ix,iy,iz) = 0.;
    B1X(fld, ix,iy,iz) =  (MRC_F3(fld_psi, 0, ix,iy+1,iz) - MRC_F3(fld_psi, 0, ix,iy,iz)) / bd2y[iy];
    B1Y(fld, ix,iy,iz) = -(MRC_F3(fld_psi, 0, ix+1,iy,iz) - MRC_F3(fld_psi, 0, ix,iy,iz)) / bd2x[ix]; 
    B1Z(fld, ix,iy,iz) = 0.;
  } mrc_fld_foreach_end;

  mrc_fld_put_as(psi, fld_psi);
  mrc_fld_put_as(fld, mhd->fld);
}

// ----------------------------------------------------------------------
// ggcm_mhd_ic_fadeev_descr

#define VAR(x) (void *)offsetof(struct ggcm_mhd_ic_fadeev, x)
static struct param ggcm_mhd_ic_fadeev_descr[] = {
  {"Bo", VAR(Bo), PARAM_FLOAT(1.0)},
  {"Boz", VAR(Boz), PARAM_FLOAT(10.0)},
  {"pert", VAR(pert), PARAM_FLOAT(0.001)},
  {"eps", VAR(eps), PARAM_FLOAT(0.3)},  
  {"lambda", VAR(lambda), PARAM_FLOAT(0.5)},  
  {"dens0", VAR(dens0), PARAM_FLOAT(5.0)},
  {},
};
#undef VAR

// ----------------------------------------------------------------------
// ggcm_mhd_ic_fadeev_ops

struct ggcm_mhd_ic_ops ggcm_mhd_ic_fadeev_ops = {
  .name        = "fadeev",
  .size        = sizeof(struct ggcm_mhd_ic_fadeev),
  .param_descr = ggcm_mhd_ic_fadeev_descr,
  .run         = ggcm_mhd_ic_fadeev_run,
};


// ======================================================================
// ggcm_mhd class "fadeev"

// ----------------------------------------------------------------------
// ggcm_mhd_fadeev_create

static void
ggcm_mhd_fadeev_create(struct ggcm_mhd *mhd)
{
  ggcm_mhd_default_box(mhd);

  ggcm_mhd_bnd_set_type(mhd->bnd, "conducting");
  mrc_domain_set_param_int(mhd->domain, "bcy", BC_NONE);	   

  /* set defaults for coord arrays */
  struct mrc_crds *crds = mrc_domain_get_crds(mhd->domain);
  mrc_crds_set_type(crds, "gaussian_2D");
  mrc_crds_set_param_int(crds, "sw", SW_2);   // 'stencil width' 
  mrc_crds_set_param_float3(crds, "l", (float[3]) {  0.0, 0.0, -1.0 });
  mrc_crds_set_param_float3(crds, "h", (float[3]) {  2.*M_PI, 2.*M_PI,  1.0 });
}

static struct ggcm_mhd_ops ggcm_mhd_fadeev_ops = {
  .name             = "fadeev",
  .create           = ggcm_mhd_fadeev_create,
};

// ======================================================================

extern struct ggcm_mhd_diag_ops ggcm_mhd_diag_c_ops;

int
main(int argc, char **argv)
{
  MPI_Init(&argc, &argv);
  libmrc_params_init(argc, argv);
  ggcm_mhd_register();

  mrc_class_register_subclass(&mrc_class_ggcm_mhd, &ggcm_mhd_fadeev_ops);  
  mrc_class_register_subclass(&mrc_class_ggcm_mhd_diag, &ggcm_mhd_diag_c_ops);

  mrc_class_register_subclass(&mrc_class_ggcm_mhd_ic, &ggcm_mhd_ic_fadeev_ops);  
 
  struct ggcm_mhd *mhd = ggcm_mhd_create(MPI_COMM_WORLD);
  ggcm_mhd_set_type(mhd, "fadeev");
  ggcm_mhd_step_set_type(mhd->step, "cweno");
  mrc_fld_set_type(mhd->fld, "float");
  ggcm_mhd_set_from_options(mhd);
  ggcm_mhd_setup(mhd);
  ggcm_mhd_view(mhd);

  // set up initial condition
  ggcm_mhd_ic_run(mhd->ic);

  ggcm_mhd_main(mhd);

  ggcm_mhd_destroy(mhd);

  MPI_Finalize();
  return 0;
}

