
#include <psc.h>

#include <mrc_params.h>

#include <math.h>

struct psc_photon_test {
  // parameters
  int n_in_cell;
  double rho_0;
  double tau_0;
  double k_L;
};

#define to_psc_photon_test(psc) mrc_to_subobj(psc, struct psc_photon_test)

#define VAR(x) (void *)offsetof(struct psc_photon_test, x)
static struct param psc_photon_test_descr[] = {
  { "n_in_cell"       , VAR(n_in_cell)         , PARAM_INT(1000)      },
  { "rho_0"           , VAR(rho_0)             , PARAM_DOUBLE(1e-6)   },
  { "tau_0"           , VAR(tau_0)             , PARAM_DOUBLE(1.5e-6) },
  { "k_L"             , VAR(k_L)               , PARAM_DOUBLE(1.)     },
  {},
};
#undef VAR

// ----------------------------------------------------------------------
// psc_photon_test_create

static void
psc_photon_test_create(struct psc *psc)
{
  psc->prm.nmax = 201;

  psc->domain.length[0] = 10 * 1e-6;
  psc->domain.length[1] = 10 * 1e-6;
  psc->domain.length[2] = 10 * 1e-6;

  psc->domain.gdims[0] = 64;
  psc->domain.gdims[1] = 64;
  psc->domain.gdims[2] = 64;
}

// ----------------------------------------------------------------------
// psc_photon_test_init_photon_np

static void
psc_photon_test_init_photon_np(struct psc *psc, double x[3],
			       struct psc_photon_np *photon_np)
{
  struct psc_photon_test *test = to_psc_photon_test(psc);

  double xc[3]; // center of blob
  double xr[3]; // distance from center, renormalized
  for (int d = 0; d < 3; d++) {
    xc[d] = .5 * ppsc->domain.length[d];
    xr[d] = x[d] * ppsc->coeff.ld - xc[0];
  };

  photon_np->n = 
    exp(-(sqr(xr[1]) + sqr(xr[2])) / sqr(test->rho_0)) *
    exp(-sqr(test->k_L * xr[0]) / sqr(test->tau_0));
  photon_np->k[0] = test->k_L;
  photon_np->sigma_k[0] = ppsc->coeff.ld / test->tau_0;
  photon_np->sigma_k[1] = ppsc->coeff.ld / test->rho_0;
  photon_np->sigma_k[2] = ppsc->coeff.ld / test->rho_0;

  // use n_in_cell particles at max of density (n == 1)
  photon_np->n_in_cell = test->n_in_cell * photon_np->n + .5;
}

// ======================================================================
// psc_photon_test_ops

struct psc_ops psc_photon_test_ops = {
  .name             = "photon_test",
  .size             = sizeof(struct psc_photon_test),
  .param_descr      = psc_photon_test_descr,
  .create           = psc_photon_test_create,
  .init_photon_np   = psc_photon_test_init_photon_np,
};

// ======================================================================
// main

int
main(int argc, char **argv)
{
  MPI_Init(&argc, &argv);
  libmrc_params_init(argc, argv);

  mrc_class_register_subclass(&mrc_class_psc, &psc_photon_test_ops);

  // psc_create() will create the psc object, create the sub-objects
  // (particle, field pusher and many others) and set the parameter defaults.
  // It will also set the psc_bubble defaults and call psc_bubble_create(),
  // which will change some of the general defaults to match this case.
  struct psc *psc = psc_create(MPI_COMM_WORLD);

  // psc_set_from_options() will override general and bubble psc parameters
  // if given on the command line. It will also call
  // psc_bubble_set_from_options()
  psc_set_from_options(psc);

  // psc_setup() will set up the various sub-objects (particle pusher, ...)
  // and set up the initial domain partition, the particles and the fields.
  // The standard implementation, used here, will set particles using
  // psc_bubble_init_npt and the fields using setup_field()
  psc_setup(psc);

  // psc_view() will just print a whole lot of info about the psc object and
  // sub-objects, in particular all the parameters.
  psc_view(psc);

  // psc_integrate() uses the standard implementation, which does the regular
  // classic PIC time integration loop
  psc_integrate(psc);

  // psc_destroy() just cleans everything up when we're done.
  psc_destroy(psc);

  libmrc_params_finalize();
  MPI_Finalize();
}
