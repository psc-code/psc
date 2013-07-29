
#include <psc.h>
#include <mrc_params.h>
#include <math.h>

// ======================================================================
// psc_em_wave
//
// This code implements a simple EM wave in vacuum

// sample command line: ~/src/psc/src/psc_em_wave --write_tfield no --pfield_step 1

#define psc_em_wave(psc) mrc_to_subobj(psc, struct psc_em_wave)

struct psc_em_wave {
  // parameters
  double k; // wave number 
};

#define VAR(x) (void *)offsetof(struct psc_em_wave, x)
static struct param psc_em_wave_descr[] = {
  { "k"             , VAR(k)                , PARAM_DOUBLE(1.)            },

  {},
};
#undef VAR

// ----------------------------------------------------------------------
// psc_em_wave_create

static void
psc_em_wave_create(struct psc *psc)
{
  psc_default_dimensionless(psc);
  psc->prm.cfl = 0.98;
  psc->prm.nmax = 100;

  psc->domain.length[0] = 1.; // no x-dependence
  psc->domain.length[1] = 1.; // no y-dependence
  psc->domain.length[2] = 2. * M_PI;

  psc->domain.gdims[0] = 1;
  psc->domain.gdims[1] = 1;
  psc->domain.gdims[2] = 32;
}

// ----------------------------------------------------------------------
// psc_em_wave_init_field

static double
psc_em_wave_init_field(struct psc *psc, double x[3], int m)
{
  struct psc_em_wave *sub = psc_em_wave(psc);
  double k = sub->k; // wave number

  switch (m) {
  case EX: return sin(k * x[2]);
  case HY: return sin(k * x[2]);
  default: return 0.;
  }
}

// ======================================================================
// psc_em_wave_ops

struct psc_ops psc_em_wave_ops = {
  .name             = "es1",
  .size             = sizeof(struct psc_em_wave),
  .param_descr      = psc_em_wave_descr,
  .create           = psc_em_wave_create,
  .init_field       = psc_em_wave_init_field,
};

// ======================================================================
// main

int
main(int argc, char **argv)
{
  return psc_main(&argc, &argv, &psc_em_wave_ops);
}
