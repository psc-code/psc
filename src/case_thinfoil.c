
#include "psc.h"
#include <mrc_params.h>

#include <math.h>
#include <string.h>
#include <stdlib.h>
#include <assert.h>

// ======================================================================
// thinfoil
//
// FIXME description

struct thinfoil {
  double Te, Ti;
  double x0, y0, z0; // location of density center in m
  double Lx, Ly, Lz; // gradient of density profile in m
  double widthx, widthy, widthz; // width of transverse / longitudinal 
                                 // density profile in m
  double rot;        // rotation angle in degrees
  double mass_ratio; // M_i / M_e
};

#define VAR(x) (void *)offsetof(struct thinfoil, x)

static struct param thinfoil_descr[] = {
  { "Te"            , VAR(Te)              , PARAM_DOUBLE(0.1)            },
  { "Ti"            , VAR(Ti)              , PARAM_DOUBLE(0.1)            },
  { "x0"            , VAR(x0)              , PARAM_DOUBLE(.01 * 1e-6)     },
  { "y0"            , VAR(y0)              , PARAM_DOUBLE(.01 * 1e-6)     },
  { "z0"            , VAR(z0)              , PARAM_DOUBLE(1.  * 1e-6)     },
  { "Lx"            , VAR(Lx)              , PARAM_DOUBLE(1.  * 1e-9)     },
  { "Ly"            , VAR(Ly)              , PARAM_DOUBLE(1.  * 1e-9)     },
  { "Lz"            , VAR(Lz)              , PARAM_DOUBLE(.5  * 1e-9)     },
  { "widthx"        , VAR(widthx)          , PARAM_DOUBLE(1.  * 1e-6)     },
  { "widthy"        , VAR(widthy)          , PARAM_DOUBLE(1.  * 1e-6)     },
  { "widthz"        , VAR(widthz)          , PARAM_DOUBLE(15. * 1e-9)     },
  { "rot"           , VAR(rot)             , PARAM_DOUBLE(0.)             },
  { "mass_ratio"    , VAR(mass_ratio)      , PARAM_DOUBLE(12.*1836.)      },
  {},
};

#undef VAR

static void
thinfoil_create(struct psc_case *Case)
{
  struct psc_pulse_flattop prm = {
    .xm = .01   * 1e-6,
    .ym = .01   * 1e-6,
    .zm = -301. * 1e-6,
    .dxm = 4.   * 1e-6,
    .dym = 4.   * 1e-6,
    .dzm = .1   * 1e-6,
    .zb  = 300. * 1e-6,
  };
  psc.pulse_p_z1 = psc_pulse_flattop_create(&prm);
}

static void
thinfoil_init_param(struct psc_case *Case)
{
  psc.prm.nmax = 10000;
  psc.prm.cpum = 25000;
  psc.prm.lw = 1. * 1e-6;
  psc.prm.i0 = 1.0e20;
  psc.prm.n0 = 1.0e29;

  psc.prm.nicell = 1000;

  psc.domain.length[0] = 0.02 * 1e-6;
  psc.domain.length[1] = 0.02 * 1e-6;
  psc.domain.length[2] = 2.0  * 1e-6;

  psc.domain.gdims[0] = 1;
  psc.domain.gdims[1] = 1;
  psc.domain.gdims[2] = 500;

  psc.domain.bnd_fld_lo[0] = 1;
  psc.domain.bnd_fld_hi[0] = 1;
  psc.domain.bnd_fld_lo[1] = 1;
  psc.domain.bnd_fld_hi[1] = 1;
  psc.domain.bnd_fld_lo[2] = 3; // time
  psc.domain.bnd_fld_hi[2] = 2; // upml
  psc.domain.bnd_part[0] = 0;
  psc.domain.bnd_part[1] = 0;
  psc.domain.bnd_part[2] = 0;
}

static void
thinfoil_init_field(struct psc_case *Case)
{
  fields_base_t *pf = &psc.pf;
  // FIXME, do we need the ghost points?
  foreach_patch(patch) {
    foreach_3d_g(patch, jx, jy, jz) {
      double dy = psc.dx[1], dz = psc.dx[2], dt = psc.dt;
      double xx = CRDX(patch, jx), yy = CRDY(patch, jy), zz = CRDZ(patch, jz);
      
      // FIXME, why this time?
      F3_BASE(pf, EY, jx,jy,jz) = psc_p_pulse_z1(xx, yy + .5*dy, zz, 0.*dt);
      F3_BASE(pf, HX, jx,jy,jz) = -psc_p_pulse_z1(xx, yy + .5*dy, zz + .5*dz, 0.*dt);
    } foreach_3d_g_end;
  }
}

static void
thinfoil_init_npt(struct psc_case *Case, int kind, double x[3], 
		  struct psc_particle_npt *npt)
{
  struct thinfoil *thinfoil = Case->ctx;

  real Te = thinfoil->Te, Ti = thinfoil->Ti;

  real ld = psc.coeff.ld;

  real x0 = thinfoil->x0 / ld;
  real y0 = thinfoil->y0 / ld;
  real z0 = thinfoil->z0 / ld;
  real Lx = thinfoil->Lx / ld;
  real Ly = thinfoil->Ly / ld;
  real Lz = thinfoil->Lz / ld;
  real widthx = thinfoil->widthx / ld;
  real widthy = thinfoil->widthy / ld;
  real widthz = thinfoil->widthz / ld;
  real rot = thinfoil->rot * 2.*M_PI / 360.;

  real xr = x[0];
  real yr = cos(rot) * (x[1]-y0) - sin(rot) * (x[2]-z0) + y0;
  real zr = cos(rot) * (x[2]-z0) + sin(rot) * (x[1]-y0) + z0;

  real argx = (fabs(xr-x0)-widthx)/Lx;
  real argy = (fabs(yr-y0)-widthy)/Ly;
  real argz = (fabs(zr-z0)-widthz)/Lz;
  if (argx > 200.0) argx = 200.0;
  if (argy > 200.0) argy = 200.0;
  if (argz > 200.0) argz = 200.0;

  real dens = 1. / ((1. + exp(argx)) * (1. + exp(argy)) * (1. + exp(argz)));

  switch (kind) {
  case 0: // electrons
    npt->q = -1.;
    npt->m = 1.;
    npt->n = dens;
    npt->T[0] = Te;
    npt->T[1] = Te;
    npt->T[2] = Te;
    break;
  case 1: // ions
    npt->q = 1.;
    npt->m = thinfoil->mass_ratio;
    npt->n = dens;
    npt->T[0] = Ti;
    npt->T[1] = Ti;
    npt->T[2] = Ti;
    break;
  default:
    assert(0);
  }
}

struct psc_case_ops psc_case_ops_thinfoil = {
  .name       = "thinfoil",
  .ctx_size   = sizeof(struct thinfoil),
  .ctx_descr  = thinfoil_descr,
  .create     = thinfoil_create,
  .init_param = thinfoil_init_param,
  .init_field = thinfoil_init_field,
  .init_npt   = thinfoil_init_npt,
};
