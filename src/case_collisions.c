
#include "psc.h"
#include "psc_case_private.h"

#include <math.h>
#include <string.h>
#include <stdlib.h>
#include <assert.h>

// ======================================================================
// collisions
//
// FIXME description

#define to_collisions(_case) ((struct collisions *)(_case)->obj.subctx)

struct collisions {
  double Te, Ti;
  double x0, y0, z0; // location of density center in m
  double Lx, Ly, Lz; // gradient of density profile in m
  double widthx, widthy, widthz; // width of transverse / longitudinal 
                                 // density profile in m
  double rot;        // rotation angle in degrees
  double mass_ratio; // M_i / M_e
};

#define VAR(x) (void *)offsetof(struct collisions, x)

static struct param collisions_descr[] = {
  { "Te"            , VAR(Te)              , PARAM_DOUBLE(0.1)             },
  { "Ti"            , VAR(Ti)              , PARAM_DOUBLE(0.1)             },
  { "x0"            , VAR(x0)              , PARAM_DOUBLE(.01 * 1e-6)      },
  { "y0"            , VAR(y0)              , PARAM_DOUBLE(.01 * 1e-6)      },
  { "z0"            , VAR(z0)              , PARAM_DOUBLE(.01 * 1e-6)      },
  { "Lx"            , VAR(Lx)              , PARAM_DOUBLE(.02 * 1e-6)      },
  { "Ly"            , VAR(Ly)              , PARAM_DOUBLE(.02 * 1e-6)      },
  { "Lz"            , VAR(Lz)              , PARAM_DOUBLE(.02 * 1e-6)      },
  { "widthx"        , VAR(widthx)          , PARAM_DOUBLE(1.  * 1e-6)      },
  { "widthy"        , VAR(widthy)          , PARAM_DOUBLE(1.  * 1e-6)      },
  { "widthz"        , VAR(widthz)          , PARAM_DOUBLE(1.  * 1e-6)      },
  { "rot"           , VAR(rot)             , PARAM_DOUBLE(0.)              },
  { "mass_ratio"    , VAR(mass_ratio)      , PARAM_DOUBLE(1836.)           },
  {},
};

#undef VAR

static void
_psc_case_collisions_set_from_options(struct _psc_case *_case)
{
  psc.prm.nmax = 10000;
  psc.prm.cpum = 25000;
  psc.prm.lw = 1. * 1e-6;
  psc.prm.i0 = 1.0e20;
  psc.prm.n0 = 0.0e35;

  psc.prm.nicell = 10000;

  psc.domain.length[0] = 0.02 * 1e-6;
  psc.domain.length[1] = 0.02 * 1e-6;
  psc.domain.length[2] = 0.02  * 1e-6;

  psc.domain.gdims[0] = 1;
  psc.domain.gdims[1] = 20;
  psc.domain.gdims[2] = 20;

  //  BND_FLD_OPEN
  //  BND_FLD_PERIODIC
  //  BND_FLD_UPML
  //  BND_FLD_TIME

  psc.domain.bnd_fld_lo[0] = BND_FLD_PERIODIC;
  psc.domain.bnd_fld_hi[0] = BND_FLD_PERIODIC;
  psc.domain.bnd_fld_lo[1] = BND_FLD_PERIODIC;
  psc.domain.bnd_fld_hi[1] = BND_FLD_PERIODIC;
  psc.domain.bnd_fld_lo[2] = BND_FLD_PERIODIC;
  psc.domain.bnd_fld_hi[2] = BND_FLD_PERIODIC;

  //  BND_PART_REFLECTING
  //  BND_PART_PERIODIC

  psc.domain.bnd_part[0] = BND_PART_PERIODIC;
  psc.domain.bnd_part[1] = BND_PART_PERIODIC;
  psc.domain.bnd_part[2] = BND_PART_PERIODIC;
}

static void
_psc_case_collisions_init_field(struct _psc_case *_case, mfields_base_t *flds)
{
  // FIXME, do we need the ghost points?
  foreach_patch(p) {
    fields_base_t *pf = &flds->f[p];
    foreach_3d_g(p, jx, jy, jz) {
      double dy = psc.dx[1], dz = psc.dx[2], dt = psc.dt;
      double xx = CRDX(p, jx), yy = CRDY(p, jy), zz = CRDZ(p, jz);
      
      // FIXME, why this time?
      F3_BASE(pf, EY, jx,jy,jz) = psc_p_pulse_z1(xx, yy + .5*dy, zz, 0.*dt);
      F3_BASE(pf, HX, jx,jy,jz) = -psc_p_pulse_z1(xx, yy + .5*dy, zz + .5*dz, 0.*dt);
    } foreach_3d_g_end;
  }
}

static void
_psc_case_collisions_init_npt(struct _psc_case *_case, int kind, double x[3], 
			      struct psc_particle_npt *npt)
{
  struct collisions *collisions = to_collisions(_case);

  real Te = collisions->Te, Ti = collisions->Ti;

  real ld = psc.coeff.ld;

  real x0 = collisions->x0 / ld;
  real y0 = collisions->y0 / ld;
  real z0 = collisions->z0 / ld;
  real Lx = collisions->Lx / ld;
  real Ly = collisions->Ly / ld;
  real Lz = collisions->Lz / ld;
  real widthx = collisions->widthx / ld;
  real widthy = collisions->widthy / ld;
  real widthz = collisions->widthz / ld;
  real rot = collisions->rot * 2.*M_PI / 360.;

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
    npt->p[1] = 0.1;
    npt->p[2] = 0.1;
    npt->T[0] = Te;
    npt->T[1] = Te;
    npt->T[2] = Te;
    break;
  case 1: // ions
    npt->q = 1.;
    npt->m = collisions->mass_ratio;
    npt->n = dens;
    npt->T[0] = Ti;
    npt->T[1] = Ti;
    npt->T[2] = Ti;
    break;
  default:
    assert(0);
  }
}

struct _psc_case_ops _psc_case_collisions_ops = {
  .name             = "collisions",
  .size             = sizeof(struct collisions),
  .param_descr      = collisions_descr,
  .set_from_options = _psc_case_collisions_set_from_options,
  .init_field       = _psc_case_collisions_init_field,
  .init_npt         = _psc_case_collisions_init_npt,
};
