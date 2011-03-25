
#include "psc.h"
#include "psc_case_private.h"

#include <math.h>
#include <string.h>
#include <stdlib.h>
#include <assert.h>

// FIXME: normalization of k etc, factors of 2 in the Gaussian (?), ...

// ----------------------------------------------------------------------
// case photon_test
//

struct psc_case_photon_test {
  // parameters
  int n_in_cell;
  double rho_0;
  double tau_0;
  double k_L;
};

#define VAR(x) (void *)offsetof(struct psc_case_photon_test, x)
static struct param psc_case_photon_test_descr[] = {
  { "n_in_cell"       , VAR(n_in_cell)         , PARAM_INT(1000)      },
  { "rho_0"           , VAR(rho_0)             , PARAM_DOUBLE(1e-6)   },
  { "tau_0"           , VAR(tau_0)             , PARAM_DOUBLE(1.5e-6) },
  { "k_L"             , VAR(k_L)               , PARAM_DOUBLE(1.)     },
  {},
};
#undef VAR

static void
psc_case_photon_test_set_from_options(struct psc_case *_case)
{
  psc.prm.nmax = 201;

  psc.domain.length[0] = 10 * 1e-6;
  psc.domain.length[1] = 10 * 1e-6;
  psc.domain.length[2] = 10 * 1e-6;

  psc.domain.gdims[0] = 64;
  psc.domain.gdims[1] = 64;
  psc.domain.gdims[2] = 64;
}

static void
psc_case_photon_test_init_photon_np(struct psc_case *_case, double x[3],
				    struct psc_photon_np *photon_np)
{
  struct psc_case_photon_test *test = mrc_to_subobj(_case, struct psc_case_photon_test);

  double xc[3]; // center of blob
  double xr[3]; // distance from center, renormalized
  for (int d = 0; d < 3; d++) {
    xc[d] = .5 * psc.domain.length[d];
    xr[d] = x[d] * psc.coeff.ld - xc[0];
  };

  photon_np->n = 
    exp(-(sqr(xr[1]) + sqr(xr[2])) / sqr(test->rho_0)) *
    exp(-sqr(test->k_L * xr[0]) / sqr(test->tau_0));
  photon_np->k[0] = test->k_L;
  photon_np->sigma_k[0] = psc.coeff.ld / test->tau_0;
  photon_np->sigma_k[1] = psc.coeff.ld / test->rho_0;
  photon_np->sigma_k[2] = psc.coeff.ld / test->rho_0;

  // use n_in_cell particles at max of density (n == 1)
  photon_np->n_in_cell = test->n_in_cell * photon_np->n + .5;
}

struct psc_case_ops psc_case_photon_test_ops = {
  .name             = "photon_test",
  .size             = sizeof(struct psc_case_photon_test),
  .param_descr      = psc_case_photon_test_descr,
  .set_from_options = psc_case_photon_test_set_from_options,
  .init_photon_np   = psc_case_photon_test_init_photon_np,
};

