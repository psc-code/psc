
#include "psc_diag.h"

static void
do_particle_energy(struct psc *psc, particles_c_t *pp, double *result)
{
  double fnqs = sqr(psc->coeff.alpha) * psc->coeff.cori / psc->coeff.eta;

  double fac = psc->dx[0] * psc->dx[1] * psc->dx[2];
  for (int n = 0; n < pp->n_part; n++) {
    particle_c_t *part = particles_c_get_one(pp, n);
      
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
}

static void
psc_diag_item_particle_energy_run(struct psc *psc, double *result)
{
  mparticles_c_t *particles = psc_mparticles_get_c(psc->particles, 0);

  psc_foreach_patch(psc, p) {
    do_particle_energy(psc, psc_mparticles_get_patch_c(particles, p), result);
  }

  psc_mparticles_put_c(particles, psc->particles); // FIXME, don't need copy-back
}

struct psc_diag_item psc_diag_item_particle_energy = {
  .run = psc_diag_item_particle_energy_run,
  .n_values = 2,
  .names = { "E_electron", "E_ion" },
};

