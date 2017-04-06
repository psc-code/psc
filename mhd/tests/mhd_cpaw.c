
#include <ggcm_mhd_private.h>
#include <ggcm_mhd_diag.h>
#include <ggcm_mhd_ic_private.h>
#include <ggcm_mhd_step.h>
#include <ggcm_mhd_defs.h>

#include <mrc_domain.h>

#include <math.h>
#include <string.h>
#include <stdio.h>
#include <stdlib.h> 
#include <assert.h>

// ======================================================================
// ggcm_mhd_ic subclass "cpaw"

struct ggcm_mhd_ic_cpaw {
  // parameters
  double rr0;
  double pp0;
  double v0;
  double B0;
  double pert;
  double k;

  // state
  double cA;
  double cs;
  double cw;
  double kx, ky;
  double amplitude_ratio;
  double alpha;
};

#define ggcm_mhd_ic_cpaw(ic) mrc_to_subobj(ic, struct ggcm_mhd_ic_cpaw)

// ----------------------------------------------------------------------
// ggcm_mhd_ic_cpaw_setup

static void
ggcm_mhd_ic_cpaw_setup(struct ggcm_mhd_ic *ic)
{
  struct ggcm_mhd_ic_cpaw *sub = ggcm_mhd_ic_cpaw(ic);
  struct ggcm_mhd *mhd = ic->mhd;
  struct mrc_crds *crds = mrc_domain_get_crds(mhd->domain);

  sub->cA = sub->B0 / sqrt(sub->rr0);
  sub->cs = sqrt(mhd->par.gamm * sub->pp0 / sub->rr0);
  const double *lo = mrc_crds_lo(crds), *hi = mrc_crds_hi(crds);
  double lambdax = hi[0] - lo[0];
  double lambday = hi[1] - lo[1];
  sub->kx = 2. * M_PI / lambdax;
  sub->ky = 2. * M_PI / lambday;

  // Initialize as strictly 1-d if domain is 1-d
  int gdims[3];
  mrc_domain_get_global_dims(ic->mhd->domain, gdims);
  if (gdims[1] == 1) {
    sub->ky = 0.;
  }

  double k = sqrt(sqr(sub->kx) + sqr(sub->ky));
  sub->alpha = atan(sub->ky / sub->kx);
  double w = mhd->par.d_i * k * fabs(sub->B0) / sub->rr0;
  sub->cw = w/2. + sqrt(sqr(sub->cA) + sqr(w/2.));
  sub->amplitude_ratio = fabs(sub->B0) / (sub->cw * sub->rr0);

  MPI_Comm comm = ggcm_mhd_ic_comm(ic);
  mpi_printf(comm, "CPAW: Alfven speed cA = %g\n", sub->cA);
  mpi_printf(comm, "CPAW: sound speed cs = %g\n", sub->cs);
  mpi_printf(comm, "CPAW: whistler speed cw = %g\n", sub->cw);
  mpi_printf(comm, "CPAW: wavelength lambda = %g\n", 2. * M_PI / k);
  mpi_printf(comm, "CPAW: wavenumber k = %g\n", k);
  mpi_printf(comm, "CPAW: rotation alpha = %g deg\n", sub->alpha * 180. / M_PI);
  mpi_printf(comm, "CPAW: amplitude ratio = %g\n", sub->amplitude_ratio);
  mpi_printf(comm, "CPAW: time for one period = %g\n", (2. * M_PI / k) / (sub->v0 + sub->cw));
}

// ----------------------------------------------------------------------
// ggcm_mhd_ic_cpaw_primitive

static double
ggcm_mhd_ic_cpaw_primitive(struct ggcm_mhd_ic *ic, int m, double crd[3])
{
  struct ggcm_mhd_ic_cpaw *sub = ggcm_mhd_ic_cpaw(ic);
  double pert = sub->pert, kx = sub->kx, ky = sub->ky, ratio = sub->amplitude_ratio;
  double alpha = sub->alpha;

  double xx = crd[0], yy = crd[1];

  switch (m) {
  case RR: return sub->rr0;
  case PP: return sub->pp0;
  case VX: return sub->v0 * cos(alpha) - pert * ratio * cos(kx * xx + ky * yy) * -sin(alpha);
  case VY: return sub->v0 * sin(alpha) - pert * ratio * cos(kx * xx + ky * yy) *  cos(alpha);
  case VZ: return                        pert * ratio * sin(kx * xx + ky * yy);
  case BX: return sub->B0 * cos(alpha) + pert * cos(kx * xx + ky * yy) * -sin(alpha);
  case BY: return sub->B0 * sin(alpha) + pert * cos(kx * xx + ky * yy) *  cos(alpha);
  case BZ: return                      - pert * sin(kx * xx + ky * yy);
  default: return 0.;
  }
}

// ----------------------------------------------------------------------
// ggcm_mhd_ic_cpaw_descr

#define VAR(x) (void *)offsetof(struct ggcm_mhd_ic_cpaw, x)
static struct param ggcm_mhd_ic_cpaw_descr[] = {
  { "rr0"         , VAR(rr0)         , PARAM_DOUBLE(1.)            },
  { "pp0"         , VAR(pp0)         , PARAM_DOUBLE(1.)            },
  { "v0"          , VAR(v0)          , PARAM_DOUBLE(0.)            },
  { "B0"          , VAR(B0)          , PARAM_DOUBLE(100.)          },
  { "pert"        , VAR(pert)        , PARAM_DOUBLE(1e-3)          },
  {},
};
#undef VAR

// ----------------------------------------------------------------------
// ggcm_mhd_ic_cpaw_ops

struct ggcm_mhd_ic_ops ggcm_mhd_ic_cpaw_ops = {
  .name        = "cpaw",
  .size        = sizeof(struct ggcm_mhd_ic_cpaw),
  .param_descr = ggcm_mhd_ic_cpaw_descr,
  .setup       = ggcm_mhd_ic_cpaw_setup,
  .primitive   = ggcm_mhd_ic_cpaw_primitive,
};


// ======================================================================
// ggcm_mhd subclass "cpaw"

// ----------------------------------------------------------------------
// ggcm_mhd_cpaw_create

static void
ggcm_mhd_cpaw_create(struct ggcm_mhd *mhd)
{
  ggcm_mhd_default_box(mhd);

  mhd->par.d_i = 35.1076;

  ggcm_mhd_step_set_type(mhd->step , "c3_double");

  /* set defaults for coord arrays */
  struct mrc_crds *crds = mrc_domain_get_crds(mhd->domain);
  mrc_crds_set_type(crds, "uniform");
  mrc_crds_set_param_int(crds, "sw", SW_2);   // 'stencil width' 
  mrc_crds_set_param_double3(crds, "l", (double[3]) {  -100., -100., -100. });
  mrc_crds_set_param_double3(crds, "h", (double[3]) {   100.,  100.,  100. });
}

// ----------------------------------------------------------------------
// ggcm_mhd_cpaw_ops

static struct ggcm_mhd_ops ggcm_mhd_cpaw_ops = {
  .name             = "cpaw",
  .create           = ggcm_mhd_cpaw_create,
};

// ======================================================================
// main

extern struct ggcm_mhd_diag_ops ggcm_mhd_diag_c_ops;

int
main(int argc, char **argv)
{
  mrc_class_register_subclass(&mrc_class_ggcm_mhd, &ggcm_mhd_cpaw_ops);  
  mrc_class_register_subclass(&mrc_class_ggcm_mhd_diag, &ggcm_mhd_diag_c_ops);
  mrc_class_register_subclass(&mrc_class_ggcm_mhd_ic, &ggcm_mhd_ic_cpaw_ops);  
 
  return ggcm_mhd_main(&argc, &argv);
}

