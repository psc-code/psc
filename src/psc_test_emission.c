
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
  psc->prm.nicell = 1;

  psc->domain.length[0] = 10 * 1e-6;
  psc->domain.length[1] = 10 * 1e-6;
  psc->domain.length[2] = 10 * 1e-6;

  psc->domain.gdims[0] = 64;
  psc->domain.gdims[1] = 64;
  psc->domain.gdims[2] = 64;
}

// ----------------------------------------------------------------------
// psc_photon_test_init_npt

static void
psc_photon_test_init_npt(struct psc *psc, int kind, double x[3],
			 struct psc_particle_npt *npt)
{
  double *dx = psc->dx;
  double xc[3];
  for (int d = 0; d < 3; d++) {
    xc[d] = psc->domain.length[d] / 2. / psc->coeff.ld;
  }

  switch (kind) {
  case 0: // electrons
    npt->q = -1.;
    npt->m = 1.;
    npt->p[0] = 2.;
    if (x[0] >= xc[0] && x[0] < xc[0] + dx[0] &&
	x[1] >= xc[1] && x[1] < xc[1] + dx[1] &&
	x[2] >= xc[2] && x[2] < xc[2] + dx[2]) {
      npt->n = 1.;
    }
    break;
  }
}

// ======================================================================
// psc_photon_test_ops

struct psc_ops psc_photon_test_ops = {
  .name             = "photon_test",
  .size             = sizeof(struct psc_photon_test),
  .param_descr      = psc_photon_test_descr,
  .create           = psc_photon_test_create,
  .init_npt         = psc_photon_test_init_npt,
};

// ======================================================================
// main

int
main(int argc, char **argv)
{
  return psc_main(&argc, &argv, &psc_photon_test_ops);
}
