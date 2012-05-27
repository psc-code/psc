
#include "psc_diag_item_private.h"
#include "psc_particles_c.h"

#include <math.h>

static void
do_particle_energy(struct psc *psc, struct psc_particles *prts_base, double *result)
{
  struct psc_particles *prts = psc_particles_get_as(prts_base, "c", 0);
  
  double fnqs = sqr(psc->coeff.alpha) * psc->coeff.cori / psc->coeff.eta;

  double fac = psc->dx[0] * psc->dx[1] * psc->dx[2];
  for (int n = 0; n < prts->n_part; n++) {
    particle_c_t *part = particles_c_get_one(prts, n);
      
    double gamma = sqrt(1.f + sqr(part->pxi) + sqr(part->pyi) + sqr(part->pzi));
    double Ekin = (gamma - 1.) * part->mni * part->wni * fnqs;
    if (part->qni < 0.) {
      result[0] += Ekin * fac;
    } else if (part->qni > 0.) {
      result[1] += Ekin * fac;
    } else {
      assert(0);
    }
  }

  psc_particles_put_as(prts, prts_base, MP_DONT_COPY);
}

static void
psc_diag_item_particle_energy_run(struct psc_diag_item *item,
				  struct psc *psc, double *result)
{
  psc_foreach_patch(psc, p) {
    do_particle_energy(psc, psc_mparticles_get_patch(psc->particles, p), result);
  }
}

// ======================================================================
// psc_diag_item_particle_energy

struct psc_diag_item_ops psc_diag_item_particle_energy_ops = {
  .name      = "particle_energy",
  .run       = psc_diag_item_particle_energy_run,
  .nr_values = 2,
  .title     = { "E_electron", "E_ion" },
};

