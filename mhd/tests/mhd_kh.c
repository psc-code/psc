
#include "ggcm_mhd_defs.h"
#include "ggcm_mhd_private.h"
#include "ggcm_mhd_crds.h"
#include "ggcm_mhd_diag.h"
#include "ggcm_mhd_ic_private.h"

#include <mrc_fld.h>
#include <mrc_fld_as_double.h>
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
  double rho0; // initial density 0 
  double rho1; // initial density 1  
  double v0; // velocity 0 
  double v1; // velocity 1
  double p0; // initial pressure
  double B0x; // initial Bx
  double B0y; // initial By
  double B0z; // initial Bz
  double B0z_harris; // initial Bz harris
  double pert; // initial pertubation amplitude
  double pert_random; // random pertubation amplitude
  double lambda; // shear layer width
  double sigma; // width of perturbation
};

// ----------------------------------------------------------------------
// ggcm_mhd_ic_kh_descr

#define VAR(x) (void *)offsetof(struct ggcm_mhd_ic_kh, x)
static struct param ggcm_mhd_ic_kh_descr[] = {
  { "rho0"          , VAR(rho0)          , PARAM_DOUBLE(2.)      },
  { "rho1"          , VAR(rho1)          , PARAM_DOUBLE(1.)      },
  { "v0"            , VAR(v0)            , PARAM_DOUBLE(.5)      },
  { "v1"            , VAR(v1)            , PARAM_DOUBLE(-.5)     },
  { "B0x"           , VAR(B0x)           , PARAM_DOUBLE(0.)      },
  { "B0y"           , VAR(B0y)           , PARAM_DOUBLE(0.)      },
  { "B0z"           , VAR(B0z)           , PARAM_DOUBLE(0.)      },
  { "B0z_harris"    , VAR(B0z_harris)    , PARAM_DOUBLE(0.)      },
  { "p0"            , VAR(p0)            , PARAM_DOUBLE(2.5)     },
  { "pert"          , VAR(pert)          , PARAM_DOUBLE(1e-2)    },
  { "pert_random"   , VAR(pert_random)   , PARAM_DOUBLE(0.)      },
  { "lambda"        , VAR(lambda)        , PARAM_DOUBLE(.05)     },
  { "sigma"         , VAR(sigma)         , PARAM_DOUBLE(.2)      },
  // FIXME, should use PARAM_SELECT
  {},
};
#undef VAR

// ----------------------------------------------------------------------
// ggcm_mhd_ic_kh0_run
//
// single mode and random perturbation, finite width shear

static void
ggcm_mhd_ic_kh0_run(struct ggcm_mhd_ic *ic)
{
  struct ggcm_mhd *mhd = ic->mhd;
  struct mrc_crds *crds = mrc_domain_get_crds(mhd->domain);  
  struct mrc_fld *fld = mrc_fld_get_as(mhd->fld, FLD_TYPE);
  struct ggcm_mhd_ic_kh *sub = mrc_to_subobj(ic, struct ggcm_mhd_ic_kh);

  int rank; MPI_Comm_rank(ggcm_mhd_ic_comm(ic), &rank);
  srandom(rank);

  for (int p = 0; p < mrc_fld_nr_patches(fld); p++) {
    mrc_fld_foreach(fld, ix,iy,iz, 1, 1) {
      double xx = MRC_MCRDX(crds, ix, p);
      double yy = MRC_MCRDY(crds, iy, p);
      
      double s = 1. + .5 * (tanh((yy - .25) / sub->lambda) - tanh((yy + .25) / sub->lambda));
      RR_(fld, ix,iy,iz, p) = sub->rho0 * s + sub->rho1 * (1. - s);
      // shear flow
      VX_(fld, ix,iy,iz, p) = sub->v0 * s + sub->v1 * (1. - s);
      // single mode perturbation
      VY_(fld, ix,iy,iz, p) = sub->pert * sin(2.*M_PI * xx) * exp(-sqr(yy) / sqr(sub->sigma));
      
      // random perturbation
      VX_(fld, ix,iy,iz, p) += sub->pert_random * (random_double() - .5);
      VY_(fld, ix,iy,iz, p) += sub->pert_random * (random_double() - .5);
      
      PP_(fld, ix,iy,iz, p) = sub->p0 
	+ .5 * sqr(sub->B0z_harris) / sqr(cosh((yy - .25) / sub->lambda))
	+ .5 * sqr(sub->B0z_harris) / sqr(cosh((yy + .25) / sub->lambda));
      BX_(fld, ix,iy,iz, p) = sub->B0x; 
      BY_(fld, ix,iy,iz, p) = sub->B0y;
      BZ_(fld, ix,iy,iz, p) = sub->B0z + sub->B0z_harris * s - sub->B0z_harris * (1. - s);
    } mrc_fld_foreach_end;
  }

  mrc_fld_put_as(fld, mhd->fld);

  ggcm_mhd_convert_from_primitive(mhd, mhd->fld);
}

// ----------------------------------------------------------------------
// ggcm_mhd_ic_kh0_ops

struct ggcm_mhd_ic_ops ggcm_mhd_ic_kh0_ops = {
  .name        = "kh",
  .size        = sizeof(struct ggcm_mhd_ic_kh),
  .param_descr = ggcm_mhd_ic_kh_descr,
  .run         = ggcm_mhd_ic_kh0_run,
};



// ----------------------------------------------------------------------
// ggcm_mhd_ic_kh1_run
//
// athena-like: random perturbations

static void
ggcm_mhd_ic_kh1_run(struct ggcm_mhd_ic *ic)
{
  struct ggcm_mhd *mhd = ic->mhd;
  struct mrc_crds *crds = mrc_domain_get_crds(mhd->domain);  
  struct mrc_fld *fld = mrc_fld_get_as(mhd->fld, FLD_TYPE);
  struct ggcm_mhd_ic_kh *sub = mrc_to_subobj(ic, struct ggcm_mhd_ic_kh);

  int rank; MPI_Comm_rank(ggcm_mhd_ic_comm(ic), &rank);
  srandom(rank);

  for (int p = 0; p < mrc_fld_nr_patches(fld); p++) {
    mrc_fld_foreach(fld, ix,iy,iz, 1, 1) {
      double yy = MRC_MCRDY(crds, iy, p);
      
      if (fabs(yy) < .25) {
	RR_(fld, ix,iy,iz, p) = sub->rho0;
	VX_(fld, ix,iy,iz, p) = sub->v0;
      } else {
	RR_(fld, ix,iy,iz, p) = sub->rho1;
	VX_(fld, ix,iy,iz, p) = sub->v1;
      }
	// random perturbation
      VX_(fld, ix,iy,iz, p) += sub->pert_random * (random_double() - .5);
      VY_(fld, ix,iy,iz, p) += sub->pert_random * (random_double() - .5);
      
      PP_(fld, ix,iy,iz, p) = sub->p0;
      
      BX_(fld, ix,iy,iz, p) = sub->B0x; 
      BY_(fld, ix,iy,iz, p) = sub->B0y;
      BZ_(fld, ix,iy,iz, p) = sub->B0z;
    } mrc_fld_foreach_end;
  }

  mrc_fld_put_as(fld, mhd->fld);

  ggcm_mhd_convert_from_primitive(mhd, mhd->fld);
}

// ----------------------------------------------------------------------
// ggcm_mhd_ic_kh1_ops

struct ggcm_mhd_ic_ops ggcm_mhd_ic_kh1_ops = {
  .name        = "kh",
  .size        = sizeof(struct ggcm_mhd_ic_kh),
  .param_descr = ggcm_mhd_ic_kh_descr,
  .run         = ggcm_mhd_ic_kh1_run,
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
  mrc_class_register_subclass(&mrc_class_ggcm_mhd_ic, &ggcm_mhd_ic_kh0_ops);  
  mrc_class_register_subclass(&mrc_class_ggcm_mhd_ic, &ggcm_mhd_ic_kh1_ops);  

  return ggcm_mhd_main(&argc, &argv);
}

