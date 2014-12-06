
#include "ggcm_mhd_defs.h"
#include "ggcm_mhd_private.h"
#include "ggcm_mhd_crds.h"
#include "ggcm_mhd_diag.h"
#include "ggcm_mhd_ic_private.h"

#include <mrc_fld.h>
#include <mrc_fld_as_float.h>
#include <mrc_domain.h>
#include <mrc_crds.h> 

#include <math.h>
#include <string.h>
#include <stdio.h>
#include <stdlib.h> 
#include <assert.h>

// ----------------------------------------------------------------------

static double
random_double()
{
  return (double) random() / RAND_MAX;
}

// ======================================================================
// ggcm_mhd_ic subclass "kh"

struct ggcm_mhd_ic_kh {
  float pert; // initial pertubation amplitude
  float pert_random; // initial random pertubation amplitude
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
  struct mrc_fld *fld = mrc_fld_get_as(mhd->fld, FLD_TYPE);
  struct ggcm_mhd_ic_kh *sub = mrc_to_subobj(ic, struct ggcm_mhd_ic_kh);
  // FIXME, the "1" no of ghosts is ugly here, and caused by the use of
  // the B1* macros which shift the index (due to staggering)...
  // FIXME, need to set all components, can't rely on things being initialized to
  // zero because of the -> primitive conversion which divides by RR :(
  double xl[3], xh[3];
  float xmid[3], L[3], r[3];
  mrc_crds_get_param_double3(crds, "l", xl);
  mrc_crds_get_param_double3(crds, "h", xh);
  for(int d = 0; d < 3; d++){
    L[d] = xh[d] - xl[d];
    xmid[d] = .5 * (xh[d] + xl[d]);
  }

  for (int p = 0; p < mrc_fld_nr_patches(fld); p++) {
    mrc_fld_foreach(fld, ix,iy,iz, 1, 1) {
      r[0] = MRC_MCRDX(crds, ix, p); 
      r[1] = MRC_MCRDY(crds, iy, p);

      float wave1 = sin(2.* sub->lambda * M_PI * (r[0] - xmid[0]) / L[0]); 
      
      if (fabs(r[1]) < .25 * L[1]) {
	RR_(fld, ix,iy,iz, p) = sub->r0;
	PP_(fld, ix,iy,iz, p) = sub->p0;
	VX_(fld, ix,iy,iz, p) = sub->v0;
	VY_(fld, ix,iy,iz, p) = sub->pert*wave1; 
      } else {
	RR_(fld, ix,iy,iz, p) = sub->r1;
	PP_(fld, ix,iy,iz, p) = sub->p0;
	VX_(fld, ix,iy,iz, p) = sub->v1;
	VY_(fld, ix,iy,iz, p) = sub->pert*wave1; 
      }
      if (sub->pert_random > 0.) {
	VX_(fld, ix,iy,iz, p) += sub->pert_random * (random_double() - .5f);
	VY_(fld, ix,iy,iz, p) += sub->pert_random * (random_double() - .5f);
      }
      BX_(fld, ix,iy,iz, p) = sub->B0; 
      BY_(fld, ix,iy,iz, p) = 0.f; 
      VZ_(fld, ix,iy,iz, p) = 0.f;
    } mrc_fld_foreach_end;
  }

  mrc_fld_put_as(fld, mhd->fld);

  ggcm_mhd_convert_from_primitive(mhd, mhd->fld);
}


// ----------------------------------------------------------------------
// ggcm_mhd_ic_kh_descr

#define VAR(x) (void *)offsetof(struct ggcm_mhd_ic_kh, x)
static struct param ggcm_mhd_ic_kh_descr[] = {
  {"pert", VAR(pert), PARAM_FLOAT(1e-2)},  
  {"pert_random", VAR(pert_random), PARAM_FLOAT(0.f)},
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
  mrc_crds_set_param_double3(crds, "l", (double[3]) {  0.0, 0.0, 0.0 });
  mrc_crds_set_param_double3(crds, "h", (double[3]) {  1.0, 1.0, 0.1 });
}

// ----------------------------------------------------------------------
// ggcm_mhd_kh_ops 

static struct ggcm_mhd_ops ggcm_mhd_kh_ops = {
  .name             = "kh",
  .create           = ggcm_mhd_kh_create,
};

// ======================================================================

extern struct ggcm_mhd_diag_ops ggcm_mhd_diag_c_ops;

int
main(int argc, char **argv)
{
  mrc_class_register_subclass(&mrc_class_ggcm_mhd, &ggcm_mhd_kh_ops);  
  mrc_class_register_subclass(&mrc_class_ggcm_mhd_diag, &ggcm_mhd_diag_c_ops);
  mrc_class_register_subclass(&mrc_class_ggcm_mhd_ic, &ggcm_mhd_ic_kh_ops);  
 
  return ggcm_mhd_main(&argc, &argv);
}

