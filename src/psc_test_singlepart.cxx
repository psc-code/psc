
#include <psc.h>
#include <psc_push_fields.h>
#include <psc_push_particles.h>
#include <psc_bnd.h>

#include <mrc_params.h>

#include <math.h>

// ======================================================================
// creates a single ion and electron, the electron is initialized to move
// at near speed of light. useful for testing pml b.c.

struct psc_test_singlepart {
  double Te, Ti;
  double loc[3]; // location of density center in m
  double electron_p[3]; // initial electron momentum
  double mass_ratio; // M_i / M_e
};

#define VAR(x) (void *)offsetof(struct psc_test_singlepart, x)
static struct param psc_test_singlepart_descr[] = {
  { "Te"            , VAR(Te)              , PARAM_DOUBLE(0.)                   },
  { "Ti"            , VAR(Ti)              , PARAM_DOUBLE(0.)                   },
  { "loc"           , VAR(loc)             , PARAM_DOUBLE3(20e-6, 20e-6, 10e-6) },
  { "electron_p"    , VAR(electron_p)      , PARAM_DOUBLE3(0., 0., 1e3)         },
  { "mass_ratio"    , VAR(mass_ratio)      , PARAM_DOUBLE(1836.)                },
  {},
};
#undef VAR

#define to_psc_test_singlepart(psc) mrc_to_subobj(psc, struct psc_test_singlepart)

// ----------------------------------------------------------------------
// psc_test_singlepart_create

static void
psc_test_singlepart_create(struct psc *psc)
{
  psc->prm.nmax = 1000;
  psc->prm.lw = 1. * 1e-6;
  psc->prm.i0 = 2.0e22;
  psc->prm.n0 = 1.0e25;

  psc->prm.nicell = 1;

  psc->domain.length[0] = 40. * 1e-6;
  psc->domain.length[1] = 40. * 1e-6;
  psc->domain.length[2] = 60. * 1e-6;

  psc->domain.gdims[0] = 1;
  psc->domain.gdims[1] = 400;
  psc->domain.gdims[2] = 600;

  psc->domain.bnd_fld_lo[0] = BND_FLD_UPML;
  psc->domain.bnd_fld_hi[0] = BND_FLD_UPML;
  psc->domain.bnd_fld_lo[1] = BND_FLD_UPML;
  psc->domain.bnd_fld_hi[1] = BND_FLD_UPML;
  psc->domain.bnd_fld_lo[2] = BND_FLD_UPML;
  psc->domain.bnd_fld_hi[2] = BND_FLD_UPML;
  psc->domain.bnd_part_lo[0] = BND_PART_REFLECTING;
  psc->domain.bnd_part_hi[0] = BND_PART_REFLECTING;
  psc->domain.bnd_part_lo[1] = BND_PART_REFLECTING;
  psc->domain.bnd_part_hi[1] = BND_PART_REFLECTING;
  psc->domain.bnd_part_lo[2] = BND_PART_REFLECTING;
  psc->domain.bnd_part_hi[2] = BND_PART_REFLECTING;
}

// ----------------------------------------------------------------------
// psc_test_singlepart_init_npt

static void
psc_test_singlepart_init_npt(struct psc *psc, int kind, double x[3],
			     struct psc_particle_npt *npt)
{
  struct psc_test_singlepart *singlepart = to_psc_test_singlepart(psc);

  real Te = singlepart->Te, Ti = singlepart->Ti;

  real loc[3];
  for (int d = 0; d < 3; d++) {
    loc[d] = singlepart->loc[d] / psc->coeff.ld;
  }
  real dens = 0.;
  double *dx = psc->patch[0].dx;
  if ((psc->domain.gdims[0] == 1 ||
       (int) (x[0] / dx[0] + .5) == (int) (loc[0] / dx[0] + .5)) && 
      (psc->domain.gdims[1] == 1 ||
       (int) (x[1] / dx[1] + .5) == (int) (loc[1] / dx[1] + .5)) && 
      (psc->domain.gdims[2] == 1 ||
       (int) (x[2] / dx[2] + .5) == (int) (loc[2] / dx[2] + .5))) {
    dens = 1.;
  }

  switch (kind) {
  case 0: // electrons
    npt->q = -1.;
    npt->m = 1.;
    npt->n = dens;
    for (int d = 0; d < 3; d++) {
      npt->p[d] = singlepart->electron_p[d];
    }
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

// ======================================================================
// psc_test_singlepart_ops

struct psc_ops psc_test_singlepart_ops = {
  .name             = "test_singlepart",
  .size             = sizeof(struct psc_test_singlepart),
  .param_descr      = psc_test_singlepart_descr,
  .create           = psc_test_singlepart_create,
  .init_npt         = psc_test_singlepart_init_npt,
};

// ======================================================================
// main

int
main(int argc, char **argv)
{
  return psc_main(&argc, &argv, &psc_test_singlepart_ops);
}
