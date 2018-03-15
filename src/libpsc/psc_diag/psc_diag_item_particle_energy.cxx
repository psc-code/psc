
#include "psc_diag_item_private.h"
#include "psc_particles_double.h"

#include <math.h>

using mparticles_t = PscMparticlesDouble;

static void
do_particle_energy(struct psc *psc, mparticles_t mprts, int p, double *result)
{
}

static void
psc_diag_item_particle_energy_run(struct psc_diag_item *item,
				  struct psc *psc, double *result)
{
  auto mprts_base = PscMparticlesBase{psc->particles};
  mparticles_t mprts = mprts_base.get_as<mparticles_t>();

  const Grid_t& grid = psc->grid();
  double fnqs = sqr(psc->coeff.alpha) * psc->coeff.cori / psc->coeff.eta;
  double fac = grid.dx[0] * grid.dx[1] * grid.dx[2];

  for (int p = 0; p < mprts->n_patches(); p++) {
    for (auto& prt : mprts[p]) {
      double gamma = sqrt(1.f + sqr(prt.pxi) + sqr(prt.pyi) + sqr(prt.pzi));
      double Ekin = (gamma - 1.) * mprts->prt_mni(prt) * mprts->prt_wni(prt) * fnqs;
      double qni = mprts->prt_qni(prt);
      if (qni < 0.) {
	result[0] += Ekin * fac;
      } else if (qni > 0.) {
	result[1] += Ekin * fac;
      } else {
	assert(0);
      }
    }
  }

  mprts.put_as(mprts_base, MP_DONT_COPY);
}

// ======================================================================
// psc_diag_item_particle_energy

struct psc_diag_item_ops_p : psc_diag_item_ops {
  psc_diag_item_ops_p() {
    name      = "particle_energy";
    run       = psc_diag_item_particle_energy_run;
    nr_values = 2;
    title[0]  = "E_electron";
    title[1]  = "E_ion";
  }
} psc_diag_item_particle_energy_ops;

