
#include "psc.h"
#include "psc_case_private.h"
#include "psc_pulse.h"
#include "psc_push_fields.h"
#include "psc_bnd_fields.h"

#include <mrc_params.h>
#include <math.h>
#include <string.h>
#include <stdlib.h>
#include <assert.h>

// ======================================================================
// thinfoil
//
// FIXME description

struct psc_case_thinfoil {
  double Te, Ti;
  double x0, y0, z0; // location of density center in m
  double Lx, Ly, Lz; // gradient of density profile in m
  double widthx, widthy, widthz; // width of transverse / longitudinal 
                                 // density profile in m
  double rot;        // rotation angle in degrees
  double mass_ratio; // M_i / M_e
};

#define VAR(x) (void *)offsetof(struct psc_case_thinfoil, x)

static struct param psc_case_thinfoil_descr[] = {
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
psc_case_thinfoil_create(struct psc_case *_case)
{
#if 0
  struct psc_pulse_flattop prm = {
    .xm = .01   * 1e-6,
    .ym = .01   * 1e-6,
    .zm = -301. * 1e-6,
    .dxm = 4.   * 1e-6,
    .dym = 4.   * 1e-6,
    .dzm = .1   * 1e-6,
    .zb  = 300. * 1e-6,
  };
#endif
  struct psc_bnd_fields *bnd_fields = psc_push_fields_get_bnd_fields(ppsc->push_fields);
  struct psc_pulse *pulse_z1 = psc_bnd_fields_get_pulse_z1(bnd_fields);
  psc_pulse_set_type(pulse_z1, "flattop");
}

static void
psc_case_thinfoil_set_from_options(struct psc_case *_case)
{
  ppsc->prm.nmax = 10000;
  ppsc->prm.lw = 1. * 1e-6;
  ppsc->prm.i0 = 1.0e20;
  ppsc->prm.n0 = 1.0e29;

  ppsc->prm.nicell = 1000;

  ppsc->domain.length[0] = 0.02 * 1e-6;
  ppsc->domain.length[1] = 0.02 * 1e-6;
  ppsc->domain.length[2] = 2.0  * 1e-6;

  ppsc->domain.gdims[0] = 1;
  ppsc->domain.gdims[1] = 1;
  ppsc->domain.gdims[2] = 500;

  ppsc->domain.bnd_fld_lo[0] = 1;
  ppsc->domain.bnd_fld_hi[0] = 1;
  ppsc->domain.bnd_fld_lo[1] = 1;
  ppsc->domain.bnd_fld_hi[1] = 1;
  ppsc->domain.bnd_fld_lo[2] = 3; // time
  ppsc->domain.bnd_fld_hi[2] = 2; // upml
  ppsc->domain.bnd_part_lo[0] = 0;
  ppsc->domain.bnd_part_hi[0] = 0;
  ppsc->domain.bnd_part_lo[1] = 0;
  ppsc->domain.bnd_part_hi[1] = 0;
  ppsc->domain.bnd_part_lo[2] = 0;
  ppsc->domain.bnd_part_hi[2] = 0;
}

static void
psc_case_thinfoil_init_npt(struct psc_case *_case, int kind, double x[3], 
			    struct psc_particle_npt *npt)
{
  struct psc_case_thinfoil *thinfoil = mrc_to_subobj(_case, struct psc_case_thinfoil);

  real Te = thinfoil->Te, Ti = thinfoil->Ti;

  real ld = ppsc->coeff.ld;

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
    npt->particles_per_cell = 100;
    break;
  default:
    assert(0);
  }
}

struct psc_case_ops psc_case_thinfoil_ops = {
  .name             = "thinfoil",
  .size             = sizeof(struct psc_case_thinfoil),
  .param_descr      = psc_case_thinfoil_descr,
  .create           = psc_case_thinfoil_create,
  .set_from_options = psc_case_thinfoil_set_from_options,
  .init_npt         = psc_case_thinfoil_init_npt,
};
