
#include "ggcm_mhd_ic_private.h"

#include "ggcm_mhd_defs.h"
#include "ggcm_mhd_private.h"
#include "ggcm_mhd_diag.h"

#include <mrc_domain.h>
#include <mrc_fld_as_float.h>
#include <mrc_crds.h>
#include <math.h>
#include <string.h>
#include <assert.h>

// ======================================================================
// ggcm_mhd_ic subclass "hydroblast"

struct ggcm_mhd_ic_hydroblast {
  float initrad; // inital radius
  float pin; // initial inside  pressure
  float pout; // initial outside pressure
  float n0; // initial density 
  const char* pdim; 
};
// ----------------------------------------------------------------------
// ggcm_mhd_ic_hydroblast_run

static void
ggcm_mhd_ic_hydroblast_run(struct ggcm_mhd_ic *ic)
{
  struct ggcm_mhd_ic_hydroblast *sub = mrc_to_subobj(ic, struct ggcm_mhd_ic_hydroblast);
  struct ggcm_mhd *mhd = ic->mhd;  
  struct mrc_fld *fld = mrc_fld_get_as(mhd->fld, FLD_TYPE);
  struct mrc_crds *crds = mrc_domain_get_crds(mhd->domain);  
  float r[3];
  for (int p = 0; p < mrc_fld_nr_patches(fld); p++) {
    mrc_fld_foreach(fld, ix, iy, iz, 1, 1) {
      r[0] = MRC_MCRD(crds, 0, ix, p);
      r[1] = MRC_MCRD(crds, 1, iy, p);
      r[2] = MRC_MCRD(crds, 2, iz, p);
      
      if((strcmp(sub->pdim, "xy") || (strcmp(sub->pdim,"yx")))== 1){
	RR_(fld, ix, iy, iz, p) = sub->n0;
	if( sqrt((r[0]*r[0]) + (r[1]*r[1])) <= sub->initrad ){
	  PP_(fld, ix, iy, iz, p) = sub->pin;
	} else{	
	  PP_(fld, ix, iy, iz, p) = sub->pout;
	}
      } else if((strcmp(sub->pdim, "yz") || (strcmp(sub->pdim,"zy"))) == 1){
	RR_(fld, ix, iy, iz, p) = sub->n0;
	if( sqrt((r[1]*r[1]) + (r[2]*r[2])) <= sub->initrad ){	
	  PP_(fld, ix, iy, iz, p) = sub->pin;
	} else{	
	  PP_(fld, ix, iy, iz, p) = sub->pout;
	}
      } else if((strcmp(sub->pdim, "xz") || (strcmp(sub->pdim,"zx"))) == 1){
	RR_(fld, ix, iy, iz, p) = sub->n0;
	if( sqrt((r[0]*r[0]) + (r[2]*r[2])) <= sub->initrad ){	
	  PP_(fld, ix, iy, iz, p) = sub->pin;
	} else {	
	  PP_(fld, ix, iy, iz, p) = sub->pout;
	}
      } else {           
	assert(0); /* unknown initial condition */
      }
    } mrc_fld_foreach_end;
  }

  mrc_fld_put_as(fld, mhd->fld);

  ggcm_mhd_convert_from_primitive(mhd, mhd->fld);
}




// ----------------------------------------------------------------------
// ggcm_mhd_ic_hydroblast_descr

#define VAR(x) (void *)offsetof(struct ggcm_mhd_ic_hydroblast, x)
static struct param ggcm_mhd_ic_hydroblast_descr[] = {
  {"initrad", VAR(initrad), PARAM_FLOAT(0.1)},
  {"pin", VAR(pin), PARAM_FLOAT(10.0)},
  {"pout", VAR(pout), PARAM_FLOAT(0.1)},
  {"n0", VAR(n0), PARAM_FLOAT(1.0)},
  {"pdim", VAR(pdim), PARAM_STRING("xy")},  
  {},
};
#undef VAR

// ----------------------------------------------------------------------
// ggcm_mhd_ic_hydroblast_ops

struct ggcm_mhd_ic_ops ggcm_mhd_ic_hydroblast_ops = {
  .name        = "hydroblast",
  .size        = sizeof(struct ggcm_mhd_ic_hydroblast),
  .param_descr = ggcm_mhd_ic_hydroblast_descr,
  .run         = ggcm_mhd_ic_hydroblast_run,
};


// ======================================================================
// ggcm_mhd class "hydroblast"

// ----------------------------------------------------------------------
// ggcm_mhd_hydroblast_create

static void
ggcm_mhd_hydroblast_create(struct ggcm_mhd *mhd)
{
  ggcm_mhd_default_box(mhd);
}

static struct ggcm_mhd_ops ggcm_mhd_hydroblast_ops = {
  .name             = "hydroblast",
  .create           = ggcm_mhd_hydroblast_create,
};

// ======================================================================
// main

extern struct ggcm_mhd_diag_ops ggcm_mhd_diag_c_ops;

int
main(int argc, char **argv)
{
  mrc_class_register_subclass(&mrc_class_ggcm_mhd, &ggcm_mhd_hydroblast_ops);  
  mrc_class_register_subclass(&mrc_class_ggcm_mhd_diag, &ggcm_mhd_diag_c_ops);
  mrc_class_register_subclass(&mrc_class_ggcm_mhd_ic, &ggcm_mhd_ic_hydroblast_ops);  
 
  return ggcm_mhd_main(&argc, &argv);
}
