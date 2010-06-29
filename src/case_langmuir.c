
#include "psc.h"

#include <math.h>
#include <string.h>
#include <stdlib.h>
#include <assert.h>

static void
langmuir_init_nvt(int kind, double x[3], double *q, double *m, double *n,
		  double v[3], double T[3])
{
  real x0 = 0.5*1.0e-6;
  real y0 = 0.5*1.0e-6;
  real z0 = 10.0*1.0e-6;
  real Lx = 0.1*1.0e-6;
  real Ly = 0.1*1.0e-6;
  real Lz = 0.1*1.0e-6;
  real widthx = 1.0*1.0e-6;
  real widthy = 1.0*1.0e-6;
  real widthz = 8.0*1.0e-6;
  real rot = 0.0;

  real ld = psc.coeff.ld;

  x0 = x0 / ld;
  y0 = y0 / ld;
  z0 = z0 / ld;
  Lx = Lx / ld;
  Ly = Ly / ld;
  Lz = Lz / ld;
  widthx = widthx / ld;
  widthy = widthy / ld;
  widthz = widthz / ld;
  rot = rot * 2.*M_PI / 360.;

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
  real dens = (.99 + .01 * cos(ld * kk * x[2])) * mask;;
  switch (kind) {
  case 0: // electrons
    *q = -1.;
    *m = 1.;
    *n = dens;
    T[2] = .05;
    break;
  case 1: // ions
    *q = 1.;
    *m = 1836.;
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
  .init_nvt   = langmuir_init_nvt,
};
