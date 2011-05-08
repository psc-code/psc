
#include "psc.h"
#include "psc_case_private.h"

#include <mrc_params.h>
#include <math.h>
#include <string.h>
#include <stdlib.h>
#include <assert.h>

// FIXME description

struct psc_case_singlepart {
  double Te, Ti;
  double x0, y0, z0; // location of density center in m
  double mass_ratio; // M_i / M_e
};

#define VAR(x) (void *)offsetof(struct psc_case_singlepart, x)

static struct param psc_case_singlepart_descr[] = {
  { "Te"            , VAR(Te)              , PARAM_DOUBLE(0.)             },
  { "Ti"            , VAR(Ti)              , PARAM_DOUBLE(0.)             },
  { "x0"            , VAR(x0)              , PARAM_DOUBLE(10.0 * 1e-6)    },
  { "y0"            , VAR(y0)              , PARAM_DOUBLE(20.0 * 1e-6)    },
  { "z0"            , VAR(z0)              , PARAM_DOUBLE(10.0 * 1e-6)    },
  { "mass_ratio"    , VAR(mass_ratio)      , PARAM_DOUBLE(1836.)          },
  {},
};

#undef VAR

static void
psc_case_singlepart_set_from_options(struct psc_case *_case)
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

  ppsc->domain.gdims[0] = 1;
  ppsc->domain.gdims[1] = 400;
  ppsc->domain.gdims[2] = 600;

  ppsc->domain.bnd_fld_lo[0] = BND_FLD_UPML;
  ppsc->domain.bnd_fld_hi[0] = BND_FLD_UPML;
  ppsc->domain.bnd_fld_lo[1] = BND_FLD_UPML;
  ppsc->domain.bnd_fld_hi[1] = BND_FLD_UPML;
  ppsc->domain.bnd_fld_lo[2] = BND_FLD_UPML;
  ppsc->domain.bnd_fld_hi[2] = BND_FLD_UPML;
  ppsc->domain.bnd_part[0] = BND_PART_REFLECTING;
  ppsc->domain.bnd_part[1] = BND_PART_REFLECTING;
  ppsc->domain.bnd_part[2] = BND_PART_REFLECTING;
}

static void
psc_case_singlepart_init_npt(struct psc_case *_case, int kind, double x[3],
			      struct psc_particle_npt *npt)
{
  struct psc_case_singlepart *singlepart = mrc_to_subobj(_case, struct psc_case_singlepart);

  real Te = singlepart->Te, Ti = singlepart->Ti;

  real ld = ppsc->coeff.ld;

  real x0[3] = { singlepart->x0 / ld,
		 singlepart->y0 / ld,
		 singlepart->z0 / ld };

  real dens = 0.0;
  if ((ppsc->domain.gdims[0] == 1 ||
       (int) (x[0] / ppsc->dx[0] + .5) == (int) (x0[0] / ppsc->dx[0] + .5)) && 
      (ppsc->domain.gdims[1] == 1 ||
       (int) (x[1] / ppsc->dx[1] + .5) == (int) (x0[1] / ppsc->dx[1] + .5)) && 
      (ppsc->domain.gdims[2] == 1 ||
       (int) (x[2] / ppsc->dx[2] + .5) == (int) (x0[2] / ppsc->dx[2] + .5))) {
    dens = 1.0;
  }

  switch (kind) {
  case 0: // electrons
    npt->q = -1.;
    npt->m = 1.;
    npt->n = dens;
    npt->p[2] = 1.0e3;
    npt->T[0] = Te;
    npt->T[1] = Te;
    npt->T[2] = Te;
    break;
  case 1: // ions
    npt->q = 1.;
    npt->m = singlepart->mass_ratio;
    npt->n = dens;
    npt->T[0] = Ti;
    npt->T[1] = Ti;
    npt->T[2] = Ti;
    break;
  default:
    assert(0);
  }
}

struct psc_case_ops psc_case_singlepart_ops = {
  .name             = "singlepart",
  .size             = sizeof(struct psc_case_singlepart),
  .param_descr      = psc_case_singlepart_descr,
  .set_from_options = psc_case_singlepart_set_from_options,
  .init_npt         = psc_case_singlepart_init_npt,
};
