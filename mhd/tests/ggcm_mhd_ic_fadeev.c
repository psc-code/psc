
#include "ggcm_mhd_ic_private.h"

#include "ggcm_mhd_defs.h"
#include "ggcm_mhd_private.h"
#include "ggcm_mhd_flds.h"
#include "ggcm_mhd_crds.h"
#include "ggcm_mhd_bnd_private.h"

#include <mrc_domain.h>
#include <mrc_crds.h>
#include <math.h>

#define ggcm_mhd_cweno(obj) mrc_to_subobj(obj, struct mhd)

// ======================================================================
// ggcm_mhd_ic subclass "fadeev"

struct ggcm_mhd_ic_fadeev {
  float mpermi;
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
  struct ggcm_mhd *gmhd = ic->mhd;
  struct mrc_f3 *f3 = ggcm_mhd_flds_get_mrc_f3(gmhd->flds_base);
  struct mrc_crds *crds = mrc_domain_get_crds(gmhd->domain);  
  struct mrc_f3 *fld_psi = mrc_domain_f3_create(gmhd->domain, SW_2);

  ggcm_mhd_bnd_set_type(gmhd->bnd, "conducting");
  //mrc_domain_set_param_int(gmhd->domain, "bcx", &bc[0]
  mrc_domain_set_param_int(gmhd->domain, "bcy", BC_NONE);	   
  mrc_f3_setup(fld_psi);

  float xl[3], xh[3], L[3], r[3];
  mrc_crds_get_xl_xh(crds, xl, xh);
  for(int i=0; i<3; i++){
    L[i] = xh[i] - xl[i];
  }
  
  float gamma = gmhd->par.gamm;
  float Bo=sub->Bo;
  float pert=sub->pert;
  float Boz=sub->Boz;
  float eps=sub->eps;
  float lam=(sub->lambda)*L[0] ;  // defines island size   
  float kk= (2.*M_PI) / lam ;      

  mrc_f3_foreach(f3, ix, iy, iz, 2, 2) {
    r[0] = .5*(MRC_CRDX(crds, ix) + MRC_CRDX(crds, ix-1));
    r[1] = .5*(MRC_CRDY(crds, iy) + MRC_CRDY(crds, iy-1));
    
    MRC_F3(fld_psi, 0, ix,iy,iz) = -(Bo / kk)*( log(cosh(kk*r[1]) + eps*cos(kk*r[0])));      
  } mrc_f3_foreach_end;

  float *bd2x = ggcm_mhd_crds_get_crd(gmhd->crds, 0, BD2);
  float *bd2y = ggcm_mhd_crds_get_crd(gmhd->crds, 1, BD2);

  mrc_f3_foreach(f3, ix, iy, iz, 1, 1) {
    // FIXME! the staggering for B is okay, but fld_psi and other stuff below needs to be
    // fixed / checked for cell-centered
    r[0] = MRC_CRD(crds, 0, ix);
    r[1] = MRC_CRD(crds, 1, iy);
    r[2] = MRC_CRD(crds, 2, iz);
    
    B1X(f3, ix,iy,iz) =  (MRC_F3(fld_psi, 0, ix,iy+1,iz) - MRC_F3(fld_psi, 0, ix,iy,iz)) / bd2y[iy];
    B1Y(f3, ix,iy,iz) = -(MRC_F3(fld_psi, 0, ix+1,iy,iz) - MRC_F3(fld_psi, 0, ix,iy,iz)) / bd2x[ix];

    MRC_F3(f3, _RR1, ix, iy, iz)  = 0.5*sqr(Bo) * (1.0-sqr(eps)) * 
      exp(2.0*kk* MRC_F3(fld_psi, 0, ix,iy,iz)/(Bo)) + 0.5*sqr(Boz) + sub->dens0;
    MRC_F3(f3, _RV1X, ix,iy,iz) = (pert) * (1.-kk*kk*r[0]*r[0]) *
      MRC_F3(f3, _RR1, ix, iy, iz) * exp(-kk*kk*r[1]*r[1])*sin(kk*r[0]*0.5);	
    MRC_F3(f3, _RV1Y, ix,iy,iz) = -(pert) * ( 0.5*kk*r[1] ) * MRC_F3(f3, _RR1, ix, iy, iz) *
      exp(-kk*kk*r[1]*r[1])*cos(kk*r[0]*0.5);            
    
    MRC_F3(f3, _UU1 , ix, iy, iz) =  MRC_F3(f3, _RR1, ix, iy, iz) / (gamma -1.f) + 	
      .5f * (sqr(MRC_F3(f3, _RV1X, ix, iy, iz)) +
	     sqr(MRC_F3(f3, _RV1Y, ix, iy, iz)) +
	     sqr(MRC_F3(f3, _RV1Z, ix, iy, iz))) / MRC_F3(f3, _RR1, ix, iy, iz) +
      .5f * (sqr(.5*(B1X(f3, ix,iy,iz) + B1X(f3, ix+1,iy,iz))) +
	     sqr(.5*(B1Y(f3, ix,iy,iz) + B1Y(f3, ix,iy+1,iz))) +
	     sqr(.5*(B1Z(f3, ix,iy,iz) + B1Z(f3, ix,iy,iz+1))));
  } mrc_f3_foreach_end;
}

// ----------------------------------------------------------------------
// ggcm_mhd_ic_fadeev_descr

#define VAR(x) (void *)offsetof(struct ggcm_mhd_ic_fadeev, x)
static struct param ggcm_mhd_ic_fadeev_descr[] = {
  {"mpermi", VAR(mpermi), PARAM_FLOAT(1.)},
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
