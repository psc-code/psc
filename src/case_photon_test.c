
#include "psc.h"
#include "psc_case_private.h"

#include <math.h>
#include <string.h>
#include <stdlib.h>
#include <assert.h>

// ----------------------------------------------------------------------
// case photon_test
//

struct psc_case_photon_test {
  // parameters
};

#define VAR(x) (void *)offsetof(struct psc_case_photon_test, x)
static struct param psc_case_photon_test_descr[] = {
  {},
};
#undef VAR

static void
psc_case_photon_test_set_from_options(struct psc_case *_case)
{
  //  struct psc_case_photon_test *msphere = mrc_to_subobj(_case, struct psc_case_photon_test);

  psc.prm.nmax = 201;

  psc.domain.length[0] = 10 * 1e-6;
  psc.domain.length[1] = 10 * 1e-6;
  psc.domain.length[2] = 10 * 1e-6;

  psc.domain.gdims[0] = 64;
  psc.domain.gdims[1] = 64;
  psc.domain.gdims[2] = 64;
}

static double
photon_test_dens(struct psc_case *_case, double x[3])
{
  // center
  const double xc[3] = { 5. * 1e-6, 5. * 1e-6, 5. * 1e-6 };
  const double w = 2 * 1e-6;

  double xr[3];
  for (int d = 0; d < 3; d++) {
    xr[d] = x[d] * psc.coeff.ld - xc[0];
  };

  double r = sqrt(sqr(xr[0]) + sqr(xr[1]) + sqr(xr[2]));

  return exp(-sqr(r / w));
}

static void
psc_case_photon_test_init_photon_np(struct psc_case *_case, double x[3],
				    struct psc_photon_np *photon_np)
{
  //  struct psc_case_photon_test *msphere = mrc_to_subobj(_case, struct psc_case_photon_test);
  double dens = photon_test_dens(_case, x);

  photon_np->n = dens;
  photon_np->p[0] = 5.;
  photon_np->n_in_cell = 10.;
}

struct psc_case_ops psc_case_photon_test_ops = {
  .name             = "photon_test",
  .size             = sizeof(struct psc_case_photon_test),
  .param_descr      = psc_case_photon_test_descr,
  .set_from_options = psc_case_photon_test_set_from_options,
  .init_photon_np   = psc_case_photon_test_init_photon_np,
};

