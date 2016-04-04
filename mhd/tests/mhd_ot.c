
#include <ggcm_mhd_private.h>
#include <ggcm_mhd_step.h>
#include <ggcm_mhd_ic_private.h>
#include <ggcm_mhd_crds_private.h>
#include <ggcm_mhd_crds_gen.h>
#include <ggcm_mhd_bnd.h>
#include <ggcm_mhd_diag.h>

#include <mrc_fld_as_double.h>
#include <mrc_domain.h>

#include <math.h>
#include <string.h>
#include <stdio.h>
#include <stdlib.h> 
#include <assert.h>

// ======================================================================
// ggcm_mhd_ic subclass "ot"

// ----------------------------------------------------------------------
// ggcm_mhd_ic_ot_primitive

static double
ggcm_mhd_ic_ot_primitive(struct ggcm_mhd_ic *ic, int m, double crd[3])
{
  double rr0 = 25. / (36.*M_PI);
  double v0 = 1.;
  double pp0 = 5. / (12.*M_PI);

  double xx = crd[0], yy = crd[1];

  switch (m) {
  case RR: return rr0;
  case PP: return pp0;
  case VX: return -v0 * sin(2. * M_PI * yy);
  case VY: return  v0 * sin(2. * M_PI * xx);
  default: return 0.;
  }
}

// ----------------------------------------------------------------------
// ggcm_mhd_ic_ot_vector_potential

static double
ggcm_mhd_ic_ot_vector_potential(struct ggcm_mhd_ic *ic, int m, double crd[3])
{
  mrc_fld_data_t B0 = 1. / sqrt(4.*M_PI);

  double xx = crd[0], yy = crd[1];

  switch (m) {
  case 2: return B0 / (4.*M_PI) * cos(4.*M_PI * xx) + B0 / (2.*M_PI) * cos(2.*M_PI * yy);
  default: return 0.;
  }
}

// ----------------------------------------------------------------------
// ggcm_mhd_ic_ot_ops

struct ggcm_mhd_ic_ops ggcm_mhd_ic_ot_ops = {
  .name             = "ot",
  .primitive        = ggcm_mhd_ic_ot_primitive,
  .vector_potential = ggcm_mhd_ic_ot_vector_potential,
};


// ======================================================================
// ggcm_mhd subclass "ot"

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
  mrc_crds_set_param_double3(crds, "l", (double[3]) {  0.0, 0.0, 0.0 });
  mrc_crds_set_param_double3(crds, "h", (double[3]) {  1.0, 1.0, 0.1 });
}

// ----------------------------------------------------------------------
// ggcm_mhd_ot_ops

static struct ggcm_mhd_ops ggcm_mhd_ot_ops = {
  .name             = "ot",
  .create           = ggcm_mhd_ot_create,
};

// ======================================================================
// main

extern struct ggcm_mhd_diag_ops ggcm_mhd_diag_c_ops;

int
main(int argc, char **argv)
{
  mrc_class_register_subclass(&mrc_class_ggcm_mhd, &ggcm_mhd_ot_ops);  
  mrc_class_register_subclass(&mrc_class_ggcm_mhd_diag, &ggcm_mhd_diag_c_ops);
  mrc_class_register_subclass(&mrc_class_ggcm_mhd_ic, &ggcm_mhd_ic_ot_ops);  
 
  return ggcm_mhd_main(&argc, &argv);
}

