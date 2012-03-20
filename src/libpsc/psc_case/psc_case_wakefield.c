
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

// FIXME description

struct psc_case_wakefield {
  double Te, Ti;
  double x0, y0, z0; // location of density center in m
  double Lx, Ly, Lz; // gradient of density profile in m
  double widthx, widthy, widthz; // width of transverse / longitudinal 
                                 // density profile in m
  double rot;        // rotation angle in degrees
  double mass_ratio; // M_i / M_e
};

#define VAR(x) (void *)offsetof(struct psc_case_wakefield, x)

static struct param psc_case_wakefield_descr[] = {
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
psc_case_wakefield_create(struct psc_case *_case)
{
  struct psc_bnd_fields *bnd_fields = psc_push_fields_get_bnd_fields(ppsc->push_fields);
  struct psc_pulse *pulse_z1 = psc_bnd_fields_get_pulse_z1(bnd_fields);
  psc_pulse_set_type(pulse_z1, "gauss");
  psc_pulse_set_param_double3(pulse_z1, "m",  (double[3]) { 10e-6, 20e-6, -2e-6 });
  psc_pulse_set_param_double3(pulse_z1, "dm", (double[3]) { 5e-6, 5e-6, 1e-6 });
  psc_pulse_set_param_double(pulse_z1, "amplitude_p", 1.);
}

static void
psc_case_wakefield_set_from_options(struct psc_case *_case)
{
  ppsc->prm.nmax = 1000;
  ppsc->prm.cpum = 20000;
  ppsc->prm.lw = 1. * 1e-6;
  ppsc->prm.i0 = 2.0e22;
  ppsc->prm.n0 = 1.0e25;

  ppsc->prm.nicell = 1;

  ppsc->domain.length[0] = 20. * 1e-6;
  ppsc->domain.length[1] = 40. * 1e-6;
  ppsc->domain.length[2] = 60. * 1e-6;

  ppsc->domain.gdims[0] = 200;
  ppsc->domain.gdims[1] = 400;
  ppsc->domain.gdims[2] = 600;

  ppsc->domain.bnd_fld_lo[0] = BND_FLD_PERIODIC;
  ppsc->domain.bnd_fld_hi[0] = BND_FLD_PERIODIC;
  ppsc->domain.bnd_fld_lo[1] = BND_FLD_PERIODIC;
  ppsc->domain.bnd_fld_hi[1] = BND_FLD_PERIODIC;
  ppsc->domain.bnd_fld_lo[2] = BND_FLD_OPEN;
  ppsc->domain.bnd_fld_hi[2] = BND_FLD_OPEN;
  ppsc->domain.bnd_part_lo[0] = BND_PART_PERIODIC;
  ppsc->domain.bnd_part_hi[0] = BND_PART_PERIODIC;
  ppsc->domain.bnd_part_lo[1] = BND_PART_PERIODIC;
  ppsc->domain.bnd_part_hi[1] = BND_PART_PERIODIC;
  ppsc->domain.bnd_part_lo[2] = BND_PART_REFLECTING;
  ppsc->domain.bnd_part_hi[2] = BND_PART_REFLECTING;
}

static void
psc_case_wakefield_init_npt(struct psc_case *_case, int kind, double x[3], 
			     struct psc_particle_npt *npt)
{
  struct psc_case_wakefield *wakefield = mrc_to_subobj(_case, struct psc_case_wakefield);

  real Te = wakefield->Te, Ti = wakefield->Ti;

  real ld = ppsc->coeff.ld;

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
    npt->T[0] = Te;
    npt->T[1] = Te;
    npt->T[2] = Te;
    break;
  case 1: // ions
    npt->q = 1.;
    npt->m = wakefield->mass_ratio;
    npt->n = dens;
    npt->T[0] = Ti;
    npt->T[1] = Ti;
    npt->T[2] = Ti;
    break;
  default:
    assert(0);
  }
}

struct psc_case_ops psc_case_wakefield_ops = {
  .name             = "wakefield",
  .size             = sizeof(struct psc_case_wakefield),
  .param_descr      = psc_case_wakefield_descr,
  .create           = psc_case_wakefield_create,
  .set_from_options = psc_case_wakefield_set_from_options,
  .init_npt         = psc_case_wakefield_init_npt,
};
