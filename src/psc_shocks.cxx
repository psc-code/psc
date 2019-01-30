
#include <psc.h>
#include <mrc_params.h>
#include <math.h>

// ======================================================================
// psc_shocks
//
// This code implements a simple EM wave in vacuum

#define psc_shocks(psc) mrc_to_subobj(psc, struct psc_shocks)

struct psc_shocks {
  // parameters
  double Tl; // left temperature
  double Tr; // right temperature
};

#define VAR(x) (void *)offsetof(struct psc_shocks, x)
static struct param psc_shocks_descr[] = {
  { "Tl"             , VAR(Tl)                , PARAM_DOUBLE(.1)            },
  { "Tr"             , VAR(Tr)                , PARAM_DOUBLE(.01)           },

  {},
};
#undef VAR

// ----------------------------------------------------------------------
// psc_shocks_create

static void
psc_shocks_create(struct psc *psc)
{
  psc_default_dimensionless(psc);
  psc->prm.cfl = 0.98;
  psc->prm.nmax = 100;

  psc->domain.length[0] = 1.; // no x-dependence
  psc->domain.length[1] = 1.; // no y-dependence
  psc->domain.length[2] = 1.;

  psc->domain.gdims[0] = 1;
  psc->domain.gdims[1] = 1;
  psc->domain.gdims[2] = 500;
}

// ----------------------------------------------------------------------
// psc_shocks_setup

static void
psc_shocks_setup(struct psc *psc)
{
  struct psc_shocks *sub = psc_shocks(psc);

  psc_setup_super(psc);

  mprintf("lambda_D = %g (l), %g (r)\n", sqrt(sub->Tl), sqrt(sub->Tr));
  double me = psc->kinds[KIND_ELECTRON].m;
  double mi = psc->kinds[KIND_ION].m;
  mprintf("d_e = %g, d_i = %g\n", sqrt(me), sqrt(mi));
}

// ----------------------------------------------------------------------
// psc_shocks_init_npt

static void
psc_shocks_init_npt(struct psc *psc, int kind, double x[3],
		    struct psc_particle_npt *npt)
{
  struct psc_shocks *sub = psc_shocks(psc);

  npt->n = 1.;
  double T;
  if (x[2] < .5 * psc->domain.length[2]) {
    T = sub->Tl;
  } else {
    T = sub->Tr;
  }

  npt->T[0] = 0.;
  npt->T[1] = 0.;
  npt->T[2] = T;
}

// ======================================================================
// psc_shocks_ops

struct psc_ops psc_shocks_ops = {
  .name             = "es1",
  .size             = sizeof(struct psc_shocks),
  .param_descr      = psc_shocks_descr,
  .create           = psc_shocks_create,
  .setup            = psc_shocks_setup,
  .init_npt         = psc_shocks_init_npt,
};

// ======================================================================
// main

int
main(int argc, char **argv)
{
  return psc_main(&argc, &argv, &psc_shocks_ops);
}

