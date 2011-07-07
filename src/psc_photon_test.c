
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
  return psc_main(&argc, &argv, &psc_photon_test_ops);
}
