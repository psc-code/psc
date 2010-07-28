
#include "psc.h"
#include "util/params.h"

#include <math.h>
#include <string.h>
#include <stdlib.h>
#include <assert.h>

// FIXME description

struct wakefield {
  double Te, Ti;
  double x0, y0, z0; // location of density center in m
  double Lx, Ly, Lz; // gradient of density profile in m
  double widthx, widthy, widthz; // width of transverse / longitudinal 
                                 // density profile in m
  double rot;        // rotation angle in degrees
  double mass_ratio; // M_i / M_e
};

#define VAR(x) (void *)offsetof(struct wakefield, x)

static struct param wakefield_descr[] = {
  { "Te"            , VAR(Te)              , PARAM_DOUBLE(0.2)    },
  { "Ti"            , VAR(Ti)              , PARAM_DOUBLE(0.)     },
  { "x0"            , VAR(x0)              , PARAM_DOUBLE(2.  * 1e-5)     },
  { "y0"            , VAR(y0)              , PARAM_DOUBLE(2.  * 1e-5)     },
  { "z0"            , VAR(z0)              , PARAM_DOUBLE(2.5 * 1e-5)     },
  { "Lx"            , VAR(Lx)              , PARAM_DOUBLE(1.  * 1e-8)     },
  { "Ly"            , VAR(Ly)              , PARAM_DOUBLE(1.  * 1e-8)     },
  { "Lz"            , VAR(Lz)              , PARAM_DOUBLE(5.  * 1e-8)     },
  { "widthx"        , VAR(widthx)          , PARAM_DOUBLE(1.  * 1e-4)     },
  { "widthy"        , VAR(widthy)          , PARAM_DOUBLE(1.  * 1e-4)     },
  { "widthz"        , VAR(widthz)          , PARAM_DOUBLE(2.  * 1e-5)     },
  { "rot"           , VAR(rot)             , PARAM_DOUBLE(0.)             },
  { "mass_ratio"    , VAR(mass_ratio)      , PARAM_DOUBLE(1836.)          },
  {},
};

#undef VAR

static void
wakefield_create(struct psc_case *Case)
{
  struct psc_p_pulse_z1_param prm = {
    .xm  = 20. * 1e-6,
    .ym  = 20. * 1e-6,
    .zm  = -2. * 1e-6,
    .dxm = 5.  * 1e-6,
    .dym = 5.  * 1e-6,
    .dzm = 1.  * 1e-6,
  };
  psc.pulse_p_z1 = psc_pulse_p_z1_short_create(&prm);
}

static void
wakefield_init_param(struct psc_case *Case)
{
  psc.prm.nmax = 1000;
  psc.prm.cpum = 20000;
  psc.prm.lw = 1. * 1e-6;
  psc.prm.i0 = 2.0e22;
  psc.prm.n0 = 1.0e25;

  psc.prm.nicell = 1;

  psc.domain.length[0] = 20. * 1e-6;
  psc.domain.length[1] = 40. * 1e-6;
  psc.domain.length[2] = 60. * 1e-6;

  psc.domain.itot[0] = 200;
  psc.domain.itot[1] = 400;
  psc.domain.itot[2] = 600;
  psc.domain.ilo[0] = 99;
  psc.domain.ilo[1] = 0;
  psc.domain.ilo[2] = 0;
  psc.domain.ihi[0] = 100;
  psc.domain.ihi[1] = 400;
  psc.domain.ihi[2] = 600;

  psc.domain.bnd_fld_lo[0] = BND_FLD_PERIODIC;
  psc.domain.bnd_fld_hi[0] = BND_FLD_PERIODIC;
  psc.domain.bnd_fld_lo[1] = BND_FLD_PERIODIC;
  psc.domain.bnd_fld_hi[1] = BND_FLD_PERIODIC;
  psc.domain.bnd_fld_lo[2] = BND_FLD_OPEN;
  psc.domain.bnd_fld_hi[2] = BND_FLD_OPEN;
  psc.domain.bnd_part[0] = BND_PART_PERIODIC;
  psc.domain.bnd_part[1] = BND_PART_PERIODIC;
  psc.domain.bnd_part[2] = BND_PART_REFLECTING;
}

static void
wakefield_init_field(struct psc_case *Case)
{
  // FIXME, do we need the ghost points?
  for (int jz = psc.ilg[2]; jz < psc.ihg[2]; jz++) {
    for (int jy = psc.ilg[1]; jy < psc.ihg[1]; jy++) {
      for (int jx = psc.ilg[0]; jx < psc.ihg[0]; jx++) {
	double dx = psc.dx[0], dy = psc.dx[1], dz = psc.dx[2], dt = psc.dt;
	double xx = jx * dx, yy = jy * dy, zz = jz * dz;

	// FIXME, why this time?
	FF3(EY, jx,jy,jz) = psc_p_pulse_z1(xx, yy + .5*dy, zz, -.5*dt);

	// FIXME, this pulse needs a - to propagate in the right direction (+z)
	FF3(BX, jx,jy,jz) = -psc_p_pulse_z1(xx, yy + .5*dy, zz + .5*dz, 0.*dt);
      }
    }
  }
}

static void
wakefield_init_npt(struct psc_case *Case, int kind, double x[3], 
		   struct psc_particle_npt *npt)
{
  struct wakefield *wakefield = Case->ctx;

  real Te = wakefield->Te, Ti = wakefield->Ti;

  real ld = psc.coeff.ld;

  real x0 = wakefield->x0 / ld;
  real y0 = wakefield->y0 / ld;
  real z0 = wakefield->z0 / ld;
  real Lx = wakefield->Lx / ld;
  real Ly = wakefield->Ly / ld;
  real Lz = wakefield->Lz / ld;
  real widthx = wakefield->widthx / ld;
  real widthy = wakefield->widthy / ld;
  real widthz = wakefield->widthz / ld;
  real rot = wakefield->rot * 2.*M_PI / 360.;

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
    npt->p[0] = 0.;
    npt->p[1] = 0.;
    npt->p[2] = 0.;
    npt->T[0] = Te;
    npt->T[1] = Te;
    npt->T[2] = Te;
    break;
  case 1: // ions
    npt->q = 1.;
    npt->m = wakefield->mass_ratio;
    npt->n = dens;
    npt->p[0] = 0.;
    npt->p[1] = 0.;
    npt->p[2] = 0.;
    npt->T[0] = Ti;
    npt->T[1] = Ti;
    npt->T[2] = Ti;
    break;
  default:
    assert(0);
  }
}

struct psc_case_ops psc_case_ops_wakefield = {
  .name       = "wakefield",
  .ctx_size   = sizeof(struct wakefield),
  .ctx_descr  = wakefield_descr,
  .create     = wakefield_create,
  .init_param = wakefield_init_param,
  .init_field = wakefield_init_field,
  .init_npt   = wakefield_init_npt,
};
