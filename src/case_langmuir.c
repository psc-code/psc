
#include "psc.h"
#include "util/params.h"

#include <math.h>
#include <string.h>
#include <stdlib.h>
#include <assert.h>

struct langmuir {
  double x0, y0, z0; // location of density center in m
  double Lx, Ly, Lz; // gradient of density profile in m
  double widthx, widthy, widthz; // width of transverse / longitudinal 
                                 // density profile in m
  double rot;        // rotation angle in degrees
  double mass_ratio; // M_i / M_e
};

#define VAR(x) (void *)offsetof(struct langmuir, x)

static struct param langmuir_descr[] = {
  { "x0"            , VAR(x0)              , PARAM_DOUBLE(.5  * 1e-6)     },
  { "y0"            , VAR(y0)              , PARAM_DOUBLE(.5  * 1e-6)     },
  { "z0"            , VAR(z0)              , PARAM_DOUBLE(10. * 1e-6)     },
  { "Lx"            , VAR(Lx)              , PARAM_DOUBLE(.1  * 1e-6)     },
  { "Ly"            , VAR(Ly)              , PARAM_DOUBLE(.1  * 1e-6)     },
  { "Lz"            , VAR(Lz)              , PARAM_DOUBLE(.1  * 1e-6)     },
  { "widthx"        , VAR(widthx)          , PARAM_DOUBLE(1.  * 1e-6)     },
  { "widthy"        , VAR(widthy)          , PARAM_DOUBLE(1.  * 1e-6)     },
  { "widthz"        , VAR(widthz)          , PARAM_DOUBLE(8.  * 1e-6)     },
  { "rot"           , VAR(rot)             , PARAM_DOUBLE(0.)             },
  { "mass_ratio"    , VAR(mass_ratio)      , PARAM_DOUBLE(1836.)          },
  {},
};

#undef VAR

static void
langmuir_create(struct psc_case *Case)
{
  struct langmuir *langmuir = Case->ctx;

  params_parse_cmdline(langmuir, langmuir_descr, "PSC Langmuir", MPI_COMM_WORLD);
  params_print(langmuir, langmuir_descr, "PSC Langmuir", MPI_COMM_WORLD);
}

static void
langmuir_init_nvt(struct psc_case *Case, 
		  int kind, double x[3], double *q, double *m, double *n,
		  double v[3], double T[3])
{
  struct langmuir *langmuir = Case->ctx;

  real ld = psc.coeff.ld;

  real x0 = langmuir->x0 / ld;
  real y0 = langmuir->y0 / ld;
  real z0 = langmuir->z0 / ld;
  real Lx = langmuir->Lx / ld;
  real Ly = langmuir->Ly / ld;
  real Lz = langmuir->Lz / ld;
  real widthx = langmuir->widthx / ld;
  real widthy = langmuir->widthy / ld;
  real widthz = langmuir->widthz / ld;
  real rot = langmuir->rot * 2.*M_PI / 360.;

  real xr = x[0];
  real yr = cos(rot) * (x[1]-y0) - sin(rot) * (x[2]-z0) + y0;
  real zr = cos(rot) * (x[2]-z0) + sin(rot) * (x[1]-y0) + z0;

  real argx = (fabs(xr-x0)-widthx)/Lx;
  real argy = (fabs(yr-y0)-widthy)/Ly;
  real argz = (fabs(zr-z0)-widthz)/Lz;
  if (argx > 200.0) argx = 200.0;
  if (argy > 200.0) argy = 200.0;
  if (argz > 200.0) argz = 200.0;

  real mask = 1. / ((1. + exp(argx)) * (1. + exp(argy)) * (1. + exp(argz)));

  real kk = 1. * M_PI * 1e6;
  real dens = (.99 + .01 * cos(ld * kk * x[2])) * mask;
  switch (kind) {
  case 0: // electrons
    *q = -1.;
    *m = 1.;
    *n = dens;
    T[2] = .05;
    break;
  case 1: // ions
    *q = 1.;
    *m = langmuir->mass_ratio;
    *n = dens;
    T[2] = 0.;
    break;
  default:
    assert(0);
  }
  v[0] = 0.;
  v[1] = 0.;
  v[2] = 0.;
  T[0] = 0.;
  T[1] = 0.;
}

struct psc_case_ops psc_case_ops_langmuir = {
  .name       = "langmuir",
  .ctx_size   = sizeof(struct langmuir),
  .ctx_descr  = langmuir_descr,
  .create     = langmuir_create,
  .init_nvt   = langmuir_init_nvt,
};
