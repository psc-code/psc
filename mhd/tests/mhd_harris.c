
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
// ggcm_mhd_ic subclass "harris"

struct ggcm_mhd_ic_harris {
  double n_0; //< peek density - n_inf
  double n_inf; //< background density
  double B_0; //< B field strength
  double cs_width; // current sheet width
  double pert; // strength of \psi perturbation (flux function)
  double pert_halfwidth; // half width of perturbation; if 0.0, set
                         // to 0.5 * (xh[1] - xl[1])
  
  // for testing: add a background field of this amount, but subtract in the perturbation,
  // so the total B field will be the same, anyway
  double B_0_bg;

  // gkeyll params
  double mass_ratio;
  double pressure_ratio;
  double d_i0;
};

#define ggcm_mhd_ic_harris(ic) mrc_to_subobj(ic, struct ggcm_mhd_ic_harris)

// ----------------------------------------------------------------------
// ggcm_mhd_ic_harris_primitive

static double
ggcm_mhd_ic_harris_primitive(struct ggcm_mhd_ic *ic, int m, double crd[3])
{
  struct ggcm_mhd_ic_harris *sub = ggcm_mhd_ic_harris(ic);
  struct ggcm_mhd *mhd = ic->mhd;
  struct mrc_crds *crds = mrc_domain_get_crds(mhd->domain);
  double pert = sub->pert, cs_width = sub->cs_width;
  double xx = crd[0], yy = crd[1];

  const double *lo = mrc_crds_lo(crds), *hi = mrc_crds_hi(crds);
  double lx = hi[0] - lo[0];
  double ly = hi[1] - lo[1];

  // FIXME, only right if cs_width == .5 (constant factor) ?
  double rr = sub->n_inf + sub->n_0 / sqr(cosh(yy/cs_width));

  switch (m) {
  case RR: return rr;
  case PP: return .5 * rr;
  // B here won't actually be used because the vector potential takes preference
  case BX: return -pert * (   M_PI/ly) * cos(2*M_PI*xx/lx) * sin(M_PI*yy/ly) + sub->B_0 * tanh(yy/cs_width)
      - sub->B_0_bg;
  case BY: return -pert * (2.*M_PI/lx) * sin(2*M_PI*xx/lx) * cos(M_PI*yy/ly);
  default: return 0.;
  }
}

// ----------------------------------------------------------------------
// ggcm_mhd_ic_harris_vector_potential

static double
ggcm_mhd_ic_harris_vector_potential(struct ggcm_mhd_ic *ic, int m, double crd[3])
{
  struct ggcm_mhd_ic_harris *sub = ggcm_mhd_ic_harris(ic);
  struct ggcm_mhd *mhd = ic->mhd;
  struct mrc_crds *crds = mrc_domain_get_crds(mhd->domain);
  double pert = sub->pert, cs_width = sub->cs_width;
  double xx = crd[0], yy = crd[1];

  const double *lo = mrc_crds_lo(crds), *hi = mrc_crds_hi(crds);
  double lx = hi[0] - lo[0];
  double ly = hi[1] - lo[1];

  switch (m) {
  case 2: return pert * cos(2*M_PI*xx/lx) * cos(M_PI*yy/ly) + sub->B_0 * cs_width * log(cosh(yy/cs_width))
      - sub->B_0_bg * yy;
  default: return 0.;
  }
}

// ----------------------------------------------------------------------
// ggcm_mhd_ic_harris_primitive_bg

static double
ggcm_mhd_ic_harris_primitive_bg(struct ggcm_mhd_ic *ic, int m, double crd[3])
{
  struct ggcm_mhd_ic_harris *sub = ggcm_mhd_ic_harris(ic);

  switch (m) {
  case BX: return sub->B_0_bg;
  default: return 0.;
  }
}

// ----------------------------------------------------------------------
// ggcm_mhd_ic_harris_vector_potential_bg

static double
ggcm_mhd_ic_harris_vector_potential_bg(struct ggcm_mhd_ic *ic, int m, double crd[3])
{
  struct ggcm_mhd_ic_harris *sub = ggcm_mhd_ic_harris(ic);
  double yy = crd[1];

  switch (m) {
  case 2: return sub->B_0_bg * yy;
  default: return 0.;
  }
}

// ----------------------------------------------------------------------
// ggcm_mhd_ic_harris_descr

#define VAR(x) (void *)offsetof(struct ggcm_mhd_ic_harris, x)
static struct param ggcm_mhd_ic_harris_descr[] = {
  { "n_0"             , VAR(n_0)             , PARAM_DOUBLE(1.0)         },
  { "n_inf"           , VAR(n_inf)           , PARAM_DOUBLE(0.2)         },
  { "B_0"             , VAR(B_0)             , PARAM_DOUBLE(1.0)         },
  { "cs_width"        , VAR(cs_width)        , PARAM_DOUBLE(0.5)         },
  { "pert"            , VAR(pert)            , PARAM_DOUBLE(0.1)         },
  { "pert_halfwidth"  , VAR(pert_halfwidth)  , PARAM_DOUBLE(6.4)         },
  { "B_0_bg"          , VAR(B_0_bg)          , PARAM_DOUBLE(0.)          },
  { "mass_ratio"      , VAR(mass_ratio)      , PARAM_DOUBLE(25.)         },
  { "pressure_ratio"  , VAR(pressure_ratio)  , PARAM_DOUBLE(1.)          },
  { "d_i0"            , VAR(d_i0)            , PARAM_DOUBLE(0.05)        },
  {},
};
#undef VAR

// ----------------------------------------------------------------------
// ggcm_mhd_ic_harris_ops

struct ggcm_mhd_ic_ops ggcm_mhd_ic_harris_ops = {
  .name                = "harris",
  .size                = sizeof(struct ggcm_mhd_ic_harris),
  .param_descr         = ggcm_mhd_ic_harris_descr,
  .primitive           = ggcm_mhd_ic_harris_primitive,
  .vector_potential    = ggcm_mhd_ic_harris_vector_potential,
  .primitive_bg        = ggcm_mhd_ic_harris_primitive_bg,
  .vector_potential_bg = ggcm_mhd_ic_harris_vector_potential_bg,
};

// ======================================================================
// ggcm_mhd_ic subclass "harris_asym"

struct ggcm_mhd_ic_harris_asym {
  // parameters
  double n0;
  double n1;
  double B0;
  double B1;
  double lambda;
  double Ti_over_Te;
  double pert; // strength of \psi perturbation (flux function)

  // state
  double nc;
  double Te, Ti;
};

#define ggcm_mhd_ic_harris_asym(ic) mrc_to_subobj(ic, struct ggcm_mhd_ic_harris_asym)

// ----------------------------------------------------------------------
// ggcm_mhd_ic_harris_asym_primitive

static double
ggcm_mhd_ic_harris_asym_primitive(struct ggcm_mhd_ic *ic, int m, double crd[3])
{
  struct ggcm_mhd_ic_harris_asym *sub = ggcm_mhd_ic_harris_asym(ic);
  double n0 = sub->n0, n1 = sub->n1, nc = sub->nc, Te = sub->Te, Ti = sub->Ti;
  double yy = crd[1];
  
  double rr = (n0 + n1) / 2. + (n0 - n1) / 2. * tanh(2.*yy) + nc / sqr(cosh(2.*yy));

  switch (m) {
  case RR: return rr;
  case PP: return (Te + Ti) * rr;
  // B is set through vector potential
  default: return 0.;
  }
}

// ----------------------------------------------------------------------
// ggcm_mhd_ic_harris_asym_vector_potential

static double
ggcm_mhd_ic_harris_asym_vector_potential(struct ggcm_mhd_ic *ic, int m, double crd[3])
{
  struct ggcm_mhd_ic_harris_asym *sub = ggcm_mhd_ic_harris_asym(ic);
  struct ggcm_mhd *mhd = ic->mhd;
  struct mrc_crds *crds = mrc_domain_get_crds(mhd->domain);
  double pert = sub->pert, B0 = sub->B0, B1 = sub->B1;
  double xx = crd[0], yy = crd[1];

  const double *lo = mrc_crds_lo(crds), *hi = mrc_crds_hi(crds);
  double lx = hi[0] - lo[0];
  double ly = hi[1] - lo[1];

  switch (m) {
  case 2: return pert * cos(2*M_PI*xx/lx) * cos(M_PI*yy/ly) + 
      (B1 - B0) / 2. * yy - (B1 + B0) / 2. * .5 * log(cosh(2.*yy));
  default: return 0.;
  }
}

// ----------------------------------------------------------------------
// ggcm_mhd_ic_harris_asym_setup

static void
ggcm_mhd_ic_harris_asym_setup(struct ggcm_mhd_ic *ic)
{
  struct ggcm_mhd_ic_harris_asym *sub = ggcm_mhd_ic_harris_asym(ic);

  sub->Te = (sqr(sub->B1) - sqr(sub->B0)) / (2. * (sub->n0 - sub->n1) * (sub->Ti_over_Te + 1.));
  sub->nc = sqr(.5 * (sub->B1 + sub->B0)) / (2. * (sub->Ti_over_Te + 1.) * sub->Te);
  sub->Ti = sub->Ti_over_Te * sub->Te;

  MPI_Comm comm = ggcm_mhd_ic_comm(ic);
  mpi_printf(comm, "i.c. harris_asym: Te = %g\n", sub->Te);
  mpi_printf(comm, "i.c. harris_asym: Ti = %g\n", sub->Ti);
  mpi_printf(comm, "i.c. harris_asym: nc = %g\n", sub->nc);
}

// ----------------------------------------------------------------------
// ggcm_mhd_ic_harris_asym_descr

#define VAR(x) (void *)offsetof(struct ggcm_mhd_ic_harris_asym, x)
static struct param ggcm_mhd_ic_harris_asym_descr[] = {
  { "n0"              , VAR(n0)              , PARAM_DOUBLE(1.)          },
  { "n1"              , VAR(n1)              , PARAM_DOUBLE(1./8.)       },
  { "B0"              , VAR(B0)              , PARAM_DOUBLE(1.)          },
  { "B1"              , VAR(B1)              , PARAM_DOUBLE(1.37)        },
  { "lambda"          , VAR(lambda)          , PARAM_DOUBLE(0.2)         },
  { "Ti_over_Te"      , VAR(Ti_over_Te)      , PARAM_DOUBLE(2.)          },
  { "pert"            , VAR(pert)            , PARAM_DOUBLE(0.1)         },
  {},
};
#undef VAR

// ----------------------------------------------------------------------
// ggcm_mhd_ic_harris_asym_ops

struct ggcm_mhd_ic_ops ggcm_mhd_ic_harris_asym_ops = {
  .name                = "harris_asym",
  .size                = sizeof(struct ggcm_mhd_ic_harris_asym),
  .param_descr         = ggcm_mhd_ic_harris_asym_descr,
  .setup               = ggcm_mhd_ic_harris_asym_setup,
  .primitive           = ggcm_mhd_ic_harris_asym_primitive,
  .vector_potential    = ggcm_mhd_ic_harris_asym_vector_potential,
};

// ======================================================================
// ggcm_mhd_ic subclass "asymharris"

struct ggcm_mhd_ic_asymharris {
  double beta_min; // min plasma beta
  double n01; //< asymtotic density
  double n02; //< asymtotic density
  double B01; //< B field strength
  double B02; //< B field strength
  double cs_width; // current sheet width
  double pert; // strength of \psi perturbation (flux function)
  double pert_halfwidth; // half width of perturbation; if 0.0, set
                         // to 0.5 * (xh[1] - xl[1])
};

// ----------------------------------------------------------------------
// ggcm_mhd_ic_asymharris_run

static void
ggcm_mhd_ic_asymharris_run(struct ggcm_mhd_ic *ic)
{
  struct ggcm_mhd_ic_asymharris *ic_asymharris = mrc_to_subobj(ic, struct ggcm_mhd_ic_asymharris);
  struct ggcm_mhd *mhd = ic->mhd;
  struct mrc_fld *fld = mrc_fld_get_as(mhd->fld, FLD_TYPE);
  struct mrc_crds *crds = mrc_domain_get_crds(mhd->domain);
  double cs_width = ic_asymharris->cs_width;

  const double *lo = mrc_crds_lo(crds), *hi = mrc_crds_hi(crds);
  double lx = hi[0] - lo[0];
  double ly = hi[1] - lo[1];
  double x0 = lo[0];
  double y0 = lo[1];
  double kx = 2.0 * M_PI / lx;
  double ky = M_PI / ly;
  double bmax = fmax(ic_asymharris->B01, ic_asymharris->B02);
  double c = (0.5 * sqr(bmax)) * (1.0 + ic_asymharris->beta_min);
  double n01 = ic_asymharris->n01;
  double n02 = ic_asymharris->n02;

  if (ic_asymharris->pert_halfwidth == 0.0) {
    ic_asymharris->pert_halfwidth = 0.5 * ly;
  }
  
  struct mrc_fld *fld_psi = mrc_domain_fld_create(fld->_domain, SW_2, "psi");
  mrc_fld_setup(fld_psi);
  
  for (int p = 0; p < mrc_fld_nr_patches(fld_psi); p++) {
    mrc_fld_foreach(fld_psi, ix,iy,iz, 2, 2) {
      double x = MRC_MCRDX(crds, ix, p);
      double y = MRC_MCRDY(crds, iy, p);

      //    A[2] = lx / (4*pi) * (1 - cos(2*kx*X)) * cos(ky*Y)
      // M3(fld_psi, 0, ix,iy,iz, p) = ic_harris->pert * ly / (4. * M_PI) * (1. - cos(2*ky*(y - y0))) * cos(kx*(x - x0));
      // taken from Birn et. al. 2001
      if (fabs(y) < ic_asymharris->pert_halfwidth) {
        M3(fld_psi, 0, ix,iy,iz, p) = (ic_asymharris->pert *
                                       cos(ky*(y - y0 - ic_asymharris->pert_halfwidth)) *
                                       cos(kx*(x - x0 - lx/2.0)));
      } else {
        M3(fld_psi, 0, ix,iy,iz, p) = 0.0;
      }
    } mrc_fld_foreach_end;
  }

  for (int p = 0; p < mrc_fld_nr_patches(fld); p++) {
    mrc_fld_foreach(fld, ix,iy,iz, 1, 2) {
      double y = MRC_MCRDY(crds, iy, p);
      double yprime = y - y0 - ly/2.0;  // y shifted to center of box
      double B0;

      if (yprime > 0.0) {
        B0 = ic_asymharris->B01;
      } else {
        B0 = ic_asymharris->B02;
      }

      // equilibrium
      BX_(fld, ix,iy,iz, p) = B0 * tanh(yprime / cs_width);
      BY_(fld, ix,iy,iz, p) = 0.0;
      BZ_(fld, ix,iy,iz, p) = 0.0;
      RR_(fld, ix,iy,iz, p) = (0.5 * (n01 + n02) +
                               0.5 * (n01 - n02) * tanh(yprime / cs_width));
      PP_(fld, ix,iy,iz, p) = c - (0.5 * sqr(BX_(fld, ix,iy,iz, p)));
      VX_(fld, ix,iy,iz, p) = 0.0;
      VY_(fld, ix,iy,iz, p) = 0.0;
      VZ_(fld, ix,iy,iz, p) = 0.0;

      // perturbation, B = -z_hat cross grad psi
      BX_(fld, ix,iy,iz, p) +=
        (M3(fld_psi, 0, ix,iy-1,iz, p) - M3(fld_psi, 0, ix,iy,iz, p)) /
        (MRC_MCRDY(crds,iy-1, p) - MRC_MCRDY(crds, iy, p));
      BY_(fld, ix,iy,iz, p) -=
        (M3(fld_psi, 0, ix-1,iy,iz, p) - M3(fld_psi, 0, ix,iy,iz, p)) /
        (MRC_MCRDX(crds,ix-1, p) - MRC_MCRDX(crds, ix, p));
    } mrc_fld_foreach_end;
  }

  mrc_fld_destroy(fld_psi);
  mrc_fld_put_as(fld, mhd->fld);

  ggcm_mhd_convert_from_primitive(mhd, mhd->fld);
}

// ----------------------------------------------------------------------
// ggcm_mhd_ic_asymharris_descr

#define VAR(x) (void *)offsetof(struct ggcm_mhd_ic_asymharris, x)
static struct param ggcm_mhd_ic_asymharris_descr[] = {
  { "beta_min"        , VAR(beta_min)        , PARAM_DOUBLE(4.0)         },
  { "n01"             , VAR(n01)             , PARAM_DOUBLE(1.0)         },
  { "n02"             , VAR(n02)             , PARAM_DOUBLE(1.0)         },
  { "B01"             , VAR(B01)             , PARAM_DOUBLE(1.0)         },
  { "B02"             , VAR(B02)             , PARAM_DOUBLE(1.0)         },
  { "cs_width"        , VAR(cs_width)        , PARAM_DOUBLE(0.5)         },
  { "pert"            , VAR(pert)            , PARAM_DOUBLE(0.1)         },
  { "pert_width"      , VAR(pert)            , PARAM_DOUBLE(6.4)         },
  {},
};
#undef VAR

// ----------------------------------------------------------------------
// ggcm_mhd_ic_asymharris_ops

struct ggcm_mhd_ic_ops ggcm_mhd_ic_asymharris_ops = {
  .name        = "asymharris",
  .size        = sizeof(struct ggcm_mhd_ic_asymharris),
  .param_descr = ggcm_mhd_ic_asymharris_descr,
  .run         = ggcm_mhd_ic_asymharris_run,
};


// ======================================================================
// ggcm_mhd subclass "harris"

// ----------------------------------------------------------------------
// ggcm_mhd_harris_create

static void
ggcm_mhd_harris_create(struct ggcm_mhd *mhd)
{
  ggcm_mhd_default_box(mhd);

  /* set defaults for coord arrays */
  struct mrc_crds *crds = mrc_domain_get_crds(mhd->domain);
  mrc_crds_set_type(crds, "uniform");
  mrc_crds_set_param_int(crds, "sw", SW_2);   // 'stencil width'
  mrc_crds_set_param_double3(crds, "l", (double[3]) {  -12.8, -6.4, 0.0 });
  mrc_crds_set_param_double3(crds, "h", (double[3]) {   12.8,  6.4, 0.1 });
  ggcm_mhd_set_param_float(mhd, "diffconstant", 0.005);

  ggcm_mhd_bnd_set_type(mhd->bnd, "conducting_y");
  mrc_domain_set_param_int(mhd->domain, "bcx", BC_PERIODIC);
  mrc_domain_set_param_int(mhd->domain, "bcy", BC_NONE);
  mrc_domain_set_param_int(mhd->domain, "bcz", BC_PERIODIC);

  ggcm_mhd_ic_set_type(mhd->ic, "harris");
}

// ----------------------------------------------------------------------
// ggcm_mhd_harris_ops

static struct ggcm_mhd_ops ggcm_mhd_harris_ops = {
  .name             = "harris",
  .create           = ggcm_mhd_harris_create,
};

// ======================================================================
// main

extern struct ggcm_mhd_diag_ops ggcm_mhd_diag_c_ops;

int
main(int argc, char **argv)
{
  mrc_class_register_subclass(&mrc_class_ggcm_mhd, &ggcm_mhd_harris_ops);
  mrc_class_register_subclass(&mrc_class_ggcm_mhd_diag, &ggcm_mhd_diag_c_ops);
  mrc_class_register_subclass(&mrc_class_ggcm_mhd_ic, &ggcm_mhd_ic_harris_ops);
  mrc_class_register_subclass(&mrc_class_ggcm_mhd_ic, &ggcm_mhd_ic_harris_asym_ops);
  mrc_class_register_subclass(&mrc_class_ggcm_mhd_ic, &ggcm_mhd_ic_asymharris_ops);

  return ggcm_mhd_main(&argc, &argv);
}
