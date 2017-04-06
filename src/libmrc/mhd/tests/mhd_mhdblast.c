
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

#define ggcm_mhd_cweno(obj) mrc_to_subobj(obj, struct mhd)

// ======================================================================
// ggcm_mhd_ic subclass "mhdblast"

struct ggcm_mhd_ic_mhdblast {
  float initrad; // inital radius
  float pin; // initial inside  pressure
  float pout; // initial outside pressure
  float n0; // initial density 
  float B0; // initial B field strentgh  
  const char* pdim; 
};
// ----------------------------------------------------------------------
// ggcm_mhd_ic_mhdblast_run

static void
ggcm_mhd_ic_mhdblast_run(struct ggcm_mhd_ic *ic)
{
  struct ggcm_mhd_ic_mhdblast *sub = mrc_to_subobj(ic, struct ggcm_mhd_ic_mhdblast);
  struct ggcm_mhd *mhd = ic->mhd;  
  struct mrc_fld *fld = mrc_fld_get_as(mhd->fld, FLD_TYPE);
  struct mrc_crds *crds = mrc_domain_get_crds(mhd->domain);  
  float r[3];
  for (int p = 0; p < mrc_fld_nr_patches(fld); p++) {
    mrc_fld_foreach(fld, ix, iy, iz, 1, 1) {
      r[0] = MRC_MCRD(crds, 0, ix, p);
      r[1] = MRC_MCRD(crds, 1, iy, p);
      r[2] = MRC_MCRD(crds, 2, iz, p);
      
      if((strcmp(sub->pdim, "xy") || (strcmp(sub->pdim,"yx"))) == 0){
	RR_(fld, ix, iy, iz, p) = sub->n0;
	BX_(fld, ix, iy, iz, p) = sub->B0/sqrt(2.0);
	BY_(fld, ix, iy, iz, p) = sub->B0/sqrt(2.0);
	if( sqrt((r[0]*r[0]) + (r[1]*r[1])) <= sub->initrad ){	
	  PP_(fld, ix, iy, iz, p) = sub->pin;
	} else{
	  PP_(fld, ix, iy, iz, p) = sub->pout;
	}
      } else if((strcmp(sub->pdim, "yz") || (strcmp(sub->pdim,"zy"))) == 0){
	RR_(fld, ix, iy, iz, p) = sub->n0;
	BX_(fld, ix, iy, iz, p) = sub->B0/sqrt(2.0);
	BY_(fld, ix, iy, iz, p) = sub->B0/sqrt(2.0);
	if( sqrt((r[0]*r[0]) + (r[1]*r[1])) <= sub->initrad ){	
	  PP_(fld, ix, iy, iz, p) = sub->pin;
	} else{
	  PP_(fld, ix, iy, iz, p) = sub->pout;
	}
      } else if((strcmp(sub->pdim, "xz") || (strcmp(sub->pdim,"zx"))) == 0){
	RR_(fld, ix, iy, iz, p) = sub->n0;
	BX_(fld, ix, iy, iz, p) = sub->B0/sqrt(2.0);
	BY_(fld, ix, iy, iz, p) = sub->B0/sqrt(2.0);
	if( sqrt((r[0]*r[0]) + (r[1]*r[1])) <= sub->initrad ){	
	  PP_(fld, ix, iy, iz, p) = sub->pin;
	} else{
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
// ggcm_mhd_ic_mhdblast_descr

#define VAR(x) (void *)offsetof(struct ggcm_mhd_ic_mhdblast, x)
static struct param ggcm_mhd_ic_mhdblast_descr[] = {
  {"pin", VAR(pin), PARAM_FLOAT(10.0)},
  {"pout", VAR(pout), PARAM_FLOAT(0.1)},
  {"n0", VAR(n0), PARAM_FLOAT(1.0)},
  {"B0", VAR(B0), PARAM_FLOAT(1.6)},
  {"pdim", VAR(pdim), PARAM_STRING("xy")},  
  {},
};
#undef VAR

// ----------------------------------------------------------------------
// ggcm_mhd_ic_mhdblast_ops

struct ggcm_mhd_ic_ops ggcm_mhd_ic_mhdblast_ops = {
  .name        = "mhdblast",
  .size        = sizeof(struct ggcm_mhd_ic_mhdblast),
  .param_descr = ggcm_mhd_ic_mhdblast_descr,
  .run         = ggcm_mhd_ic_mhdblast_run,
};


// ======================================================================
// ggcm_mhd class "mhdblast"

// ----------------------------------------------------------------------
// ggcm_mhd_mhdblast_create

static void
ggcm_mhd_mhdblast_create(struct ggcm_mhd *mhd)
{
  ggcm_mhd_default_box(mhd);
}

static struct ggcm_mhd_ops ggcm_mhd_mhdblast_ops = {
  .name             = "mhdblast",
  .create           = ggcm_mhd_mhdblast_create,
};

// ======================================================================
// main

extern struct ggcm_mhd_diag_ops ggcm_mhd_diag_c_ops;

int
main(int argc, char **argv)
{
  mrc_class_register_subclass(&mrc_class_ggcm_mhd, &ggcm_mhd_mhdblast_ops);  
  mrc_class_register_subclass(&mrc_class_ggcm_mhd_diag, &ggcm_mhd_diag_c_ops);
  mrc_class_register_subclass(&mrc_class_ggcm_mhd_ic, &ggcm_mhd_ic_mhdblast_ops);  
 
  return ggcm_mhd_main(&argc, &argv);
}
