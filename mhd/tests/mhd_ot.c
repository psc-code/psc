
#include <ggcm_mhd_private.h>
#include <ggcm_mhd_step.h>
#include <ggcm_mhd_ic_private.h>
#include <ggcm_mhd_crds_private.h>
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

struct ggcm_mhd_ic_ot {
  // params
  double rr0;
  double v0;
  double pp0;
  double B0;

  // gkeyll params
  double mass_ratio;
  double pressure_ratio;
  double d_i0;

  // state
  double kx;
  double ky;
};

#define ggcm_mhd_ic_ot(ic) mrc_to_subobj(ic, struct ggcm_mhd_ic_ot)

// ----------------------------------------------------------------------
// ggcm_mhd_ic_ot_setup

static void
ggcm_mhd_ic_ot_setup(struct ggcm_mhd_ic *ic)
{
  struct ggcm_mhd_ic_ot *sub = ggcm_mhd_ic_ot(ic);
  struct mrc_crds *crds = mrc_domain_get_crds(ic->mhd->domain);
  const double *lo = mrc_crds_lo(crds), *hi = mrc_crds_hi(crds);
  
  sub->kx = 2. * M_PI / (hi[0] - lo[0]);
  sub->ky = 2. * M_PI / (hi[1] - lo[1]);
}

// ----------------------------------------------------------------------
// ggcm_mhd_ic_ot_primitive

static double
ggcm_mhd_ic_ot_primitive(struct ggcm_mhd_ic *ic, int m, double crd[3])
{
  struct ggcm_mhd_ic_ot *sub = ggcm_mhd_ic_ot(ic);

  double rr0 = sub->rr0, pp0 = sub->pp0, v0 = sub->v0;
  double B0 = sub->B0, kx = sub->kx, ky = sub->ky;
  double xx = crd[0], yy = crd[1];

  switch (m) {
  case RR: return rr0;
  case PP: return pp0;
  case VX: return -v0 * sin(ky * yy);
  case VY: return  v0 * sin(kx * xx);
  // B here won't actually be used because the vector potential takes preference
  case BX: return -B0 * sin(ky * yy);
  case BY: return  B0 * sin(2. * kx * xx);
  default: return 0.;
  }
}

// ----------------------------------------------------------------------
// ggcm_mhd_ic_ot_vector_potential

static double
ggcm_mhd_ic_ot_vector_potential(struct ggcm_mhd_ic *ic, int m, double crd[3])
{
  struct ggcm_mhd_ic_ot *sub = ggcm_mhd_ic_ot(ic);

  double B0 = sub->B0, kx = sub->kx, ky = sub->ky;
  double xx = crd[0], yy = crd[1];

  switch (m) {
  case 2: return B0 / (2. * kx) * cos(2. * kx * xx) + B0 / ky * cos(ky * yy);
  default: return 0.;
  }
}

// ----------------------------------------------------------------------
// ggcm_mhd_ic_ot_descr

#define VAR(x) (void *)offsetof(struct ggcm_mhd_ic_ot, x)
static struct param ggcm_mhd_ic_ot_descr[] = {
  { "rr0"           , VAR(rr0)           , PARAM_DOUBLE(25. / (36. * M_PI))   },
  { "v0"            , VAR(v0)            , PARAM_DOUBLE(1.)                   },
  { "pp0"           , VAR(pp0)           , PARAM_DOUBLE(5. / (12. * M_PI))    },
  { "B0"            , VAR(B0)            , PARAM_DOUBLE(0.28209479177387814)  }, // 1. / sqrt(4. * M_PI)
  { "mass_ratio"    , VAR(mass_ratio)    , PARAM_DOUBLE(25.)                  },
  { "pressure_ratio", VAR(pressure_ratio), PARAM_DOUBLE(1.)                   },
  { "d_i0"          , VAR(d_i0)          , PARAM_DOUBLE(0.05)                 },
  {},
};
#undef VAR

// ----------------------------------------------------------------------
// ggcm_mhd_ic_ot_ops

struct ggcm_mhd_ic_ops ggcm_mhd_ic_ot_ops = {
  .name             = "ot",
  .size             = sizeof(struct ggcm_mhd_ic_ot),
  .param_descr      = ggcm_mhd_ic_ot_descr,
  .setup            = ggcm_mhd_ic_ot_setup,
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

  // default mesh size
  mrc_domain_set_param_int3(mhd->domain, "m", (int[3]) { 64, 64, 1 });

  // default domain size
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

