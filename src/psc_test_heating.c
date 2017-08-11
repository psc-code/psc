
#include <psc.h>
#include <psc_push_particles.h>
#include <psc_push_fields.h>
#include <psc_sort.h>
#include <psc_balance.h>
#include <psc_particles_as_single.h>

#include <mrc_params.h>
#include <mrc_profile.h>

#include <math.h>
#include <time.h>

// ======================================================================

struct psc_test_heating {
  // parameters
  double mi_over_me;
  double wpe_over_wce;
  double Ti_over_Te;
  double beta;

  // calculated from the above
  double B0;
};

#define psc_test_heating(psc) mrc_to_subobj(psc, struct psc_test_heating)

#define VAR(x) (void *)offsetof(struct psc_test_heating, x)
static struct param psc_test_heating_descr[] = {
  { "mi_over_me"    , VAR(mi_over_me)      , PARAM_DOUBLE(5.)            },
  { "wpe_over_wce"  , VAR(wpe_over_wce)    , PARAM_DOUBLE(2.)            },
  { "Ti_over_Te"    , VAR(Ti_over_Te)      , PARAM_DOUBLE(1.)            },
  { "beta"          , VAR(beta)            , PARAM_DOUBLE(.5)            },
  {},
};
#undef VAR

// ----------------------------------------------------------------------
// psc_test_heating_create

static void
psc_test_heating_create(struct psc *psc)
{
  psc_default_dimensionless(psc);

  psc->prm.nmax = 16000;
  psc->prm.nicell = 50;
  psc->prm.cfl = 0.98;

  psc->domain.length[0] = 1.; // no x-dependence
  psc->domain.length[1] = 30.;
  psc->domain.length[2] = 30.;

  psc->domain.gdims[0] = 1;
  psc->domain.gdims[1] = 160;
  psc->domain.gdims[2] = 160;

  psc->domain.bnd_fld_lo[0] = BND_FLD_PERIODIC;
  psc->domain.bnd_fld_hi[0] = BND_FLD_PERIODIC;
  psc->domain.bnd_fld_lo[1] = BND_FLD_PERIODIC;
  psc->domain.bnd_fld_hi[1] = BND_FLD_PERIODIC;
  psc->domain.bnd_fld_lo[2] = BND_FLD_PERIODIC;
  psc->domain.bnd_fld_hi[2] = BND_FLD_PERIODIC;

  psc->domain.bnd_part_lo[0] = BND_PART_PERIODIC;
  psc->domain.bnd_part_hi[0] = BND_PART_PERIODIC;
  psc->domain.bnd_part_lo[1] = BND_PART_PERIODIC;
  psc->domain.bnd_part_hi[1] = BND_PART_PERIODIC;
  psc->domain.bnd_part_lo[2] = BND_PART_PERIODIC;
  psc->domain.bnd_part_hi[2] = BND_PART_PERIODIC;
}

// ----------------------------------------------------------------------
// psc_test_heating_setup
//
// the parameters are now set, calculate quantities to initialize fields,
// particles

static void
psc_test_heating_setup(struct psc *psc)
{
  struct psc_test_heating *sub = psc_test_heating(psc);

  double me = 1.;
  double mi = me * sub->mi_over_me;
  double B0 = sqrt(me) / (sub->wpe_over_wce);
  double Te = sub->beta * (1. / (1. + sub->Ti_over_Te)) * sqr(B0) / 2.;
  double Ti = sub->beta * (1. / (1. + 1./sub->Ti_over_Te)) * sqr(B0) / 2.;
  mpi_printf(MPI_COMM_WORLD, "psc/test_heating: B0 %g\n", B0);
  mpi_printf(MPI_COMM_WORLD, "psc/test_heating: Te %g Ti %g\n", Te, Ti);
  mpi_printf(MPI_COMM_WORLD, "psc/test_heating: ve %g vi %g\n", sqrt(Te/me), sqrt(Ti/mi));
  mpi_printf(MPI_COMM_WORLD, "psc/test_heating: lambda_De %g\n", sqrt(Te));

  sub->B0 = B0;

  psc->kinds[KIND_ELECTRON].m = me;
  psc->kinds[KIND_ELECTRON].T = Te;
  psc->kinds[KIND_ION].m = mi;
  psc->kinds[KIND_ION].T = Ti;

  psc_setup_super(psc);
}

// ----------------------------------------------------------------------
// psc_test_heating_init_field

static double
psc_test_heating_init_field(struct psc *psc, double x[3], int m)
{
  struct psc_test_heating *sub = psc_test_heating(psc);

  switch (m) {
  case HZ: return sub->B0;
  default: return 0.;
  }
}

// ----------------------------------------------------------------------
// psc_test_heating_init_npt

static void
psc_test_heating_init_npt(struct psc *psc, int kind, double x[3],
			  struct psc_particle_npt *npt)
{
  npt->n = 1;
}

// ======================================================================
// psc_test_heating_ops

struct psc_ops psc_test_heating_ops = {
  .name             = "test_heating",
  .size             = sizeof(struct psc_test_heating),
  .param_descr      = psc_test_heating_descr,
  .create           = psc_test_heating_create,
  .setup            = psc_test_heating_setup,
  .init_field       = psc_test_heating_init_field,
  .init_npt         = psc_test_heating_init_npt,
};

// ======================================================================
// main

int
main(int argc, char **argv)
{
  return psc_main(&argc, &argv, &psc_test_heating_ops);
}
