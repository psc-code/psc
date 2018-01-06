
#include <psc.h>
#include <psc_push_particles.h>
#include <psc_push_fields.h>
#include <psc_collision.h>
#include <psc_balance.h>

#include <mrc_params.h>

#include <math.h>

struct psc_test_open {
  double B0;
  double mi_over_me;
  double Te, Ti;
  double vye, vze;
  double vyi, vzi;
};

#define to_psc_test_open(psc) mrc_to_subobj(psc, struct psc_test_open)

#define VAR(x) (void *)offsetof(struct psc_test_open, x)
static struct param psc_test_open_descr[] = {
  { "B0"            , VAR(B0)              , PARAM_DOUBLE(0.)     },
  { "mi_over_me"    , VAR(mi_over_me)      , PARAM_DOUBLE(25.)    },
  { "Te"            , VAR(Te)              , PARAM_DOUBLE(.1)     },
  { "Ti"            , VAR(Ti)              , PARAM_DOUBLE(.1)     },
  { "vye"           , VAR(vye)             , PARAM_DOUBLE(0.)     },
  { "vze"           , VAR(vze)             , PARAM_DOUBLE(.1)     },
  { "vyi"           , VAR(vyi)             , PARAM_DOUBLE(0.)     },
  { "vzi"           , VAR(vzi)             , PARAM_DOUBLE(.1)     },
  {},
};
#undef VAR

// ----------------------------------------------------------------------
// psc_test_open_create

static void
psc_test_open_create(struct psc *psc)
{
  psc_default_dimensionless(psc);

  psc->prm.nmax = 16000;
  psc->prm.nicell = 50;
  psc->prm.cfl = 0.98;

  // will be set to actual values in psc_test_open_setup()
  psc->domain.length[0] = 1.; // no x dependence 
  psc->domain.length[1] = 6.4;
  psc->domain.length[2] = 12.8;

  psc->domain.gdims[0] = 1;
  psc->domain.gdims[1] = 64;
  psc->domain.gdims[2] = 128;

  psc->domain.bnd_fld_lo[0] = BND_FLD_PERIODIC;
  psc->domain.bnd_fld_hi[0] = BND_FLD_PERIODIC;
  psc->domain.bnd_fld_lo[1] = BND_FLD_OPEN;
  psc->domain.bnd_fld_hi[1] = BND_FLD_OPEN;
  psc->domain.bnd_fld_lo[2] = BND_FLD_OPEN;
  psc->domain.bnd_fld_hi[2] = BND_FLD_OPEN;
 
  psc->domain.bnd_part_lo[0] = BND_PART_PERIODIC;
  psc->domain.bnd_part_hi[0] = BND_PART_PERIODIC;
  psc->domain.bnd_part_lo[1] = BND_PART_OPEN;
  psc->domain.bnd_part_hi[1] = BND_PART_OPEN;
  psc->domain.bnd_part_lo[2] = BND_PART_OPEN;
  psc->domain.bnd_part_hi[2] = BND_PART_OPEN;

  // FIXME: can only use 1st order pushers with current conducting wall b.c.
  psc_push_particles_set_type(psc->push_particles, "1vb");
}

// ----------------------------------------------------------------------
// psc_test_open_setup

static void
psc_test_open_setup(struct psc *psc)
{
  struct psc_test_open *sub = to_psc_test_open(psc);

  psc->kinds[KIND_ELECTRON].m = 1.;
  psc->kinds[KIND_ION     ].m = sub->mi_over_me;
  psc->kinds[KIND_ELECTRON].T = sub->Te;
  psc->kinds[KIND_ION     ].T = sub->Ti;

  // initializes fields, particles, etc.
  psc_setup_super(psc);
}

// ----------------------------------------------------------------------
// psc_test_open_init_field

static double
psc_test_open_init_field(struct psc *psc, double x[3], int m)
{
  //  struct psc_test_open *sub = to_psc_test_open(psc);

  //double B0 = sub->B0;

  return 0.;
}

// ----------------------------------------------------------------------
// psc_test_open_init_npt
//
// jx = n e (vi - ve) = 1/cosh^2 (2 (TTi + TTe) / BB / LLL

static void
psc_test_open_init_npt(struct psc *psc, int pop, double x[3],
		       struct psc_particle_npt *npt)
{
  struct psc_test_open *sub = to_psc_test_open(psc);

  npt->n = 1.;
  switch (pop) {
  case KIND_ELECTRON:
    npt->q = -1.;
    npt->m = 1.;
    npt->p[1] = sub->vye;
    npt->p[2] = sub->vze;
    break;
  case KIND_ION:
    npt->q = 1.;
    npt->m = sub->mi_over_me;
    npt->p[1] = sub->vyi;
    npt->p[2] = sub->vzi;
    break;
  default:
    assert(0);
  }
}

// ======================================================================
// psc_test_open_ops

struct psc_ops psc_test_open_ops = {
  .name             = "test_open",
  .size             = sizeof(struct psc_test_open),
  .param_descr      = psc_test_open_descr,
  .create           = psc_test_open_create,
  .setup            = psc_test_open_setup,
  .init_field       = psc_test_open_init_field,
  .init_npt         = psc_test_open_init_npt,
};

// ======================================================================
// main

int
main(int argc, char **argv)
{
  return psc_main(&argc, &argv, &psc_test_open_ops);
}
