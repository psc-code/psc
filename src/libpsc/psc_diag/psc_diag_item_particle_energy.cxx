
#include "psc_diag_item_private.h"
#include "psc_particles_as_double.h"

#include <math.h>

static void
do_particle_energy(struct psc *psc, mparticles_t mprts, int p, double *result)
{
  particle_range_t prts = mprts.range(p);
  double fnqs = sqr(psc->coeff.alpha) * psc->coeff.cori / psc->coeff.eta;

  struct psc_patch *patch = &psc->patch[p];
  double fac = patch->dx[0] * patch->dx[1] * patch->dx[2];
  for (int n = 0; n < particle_range_size(prts); n++) {
    particle_t *part = particle_iter_at(prts.begin, n);
      
    double gamma = sqrt(1.f + sqr(part->pxi) + sqr(part->pyi) + sqr(part->pzi));
    double Ekin = (gamma - 1.) * particle_mni(part) * particle_wni(part) * fnqs;
    if (particle_qni(part) < 0.) {
      result[0] += Ekin * fac;
    } else if (particle_qni(part) > 0.) {
      result[1] += Ekin * fac;
    } else {
      assert(0);
    }
  }
}

static void
psc_diag_item_particle_energy_run(struct psc_diag_item *item,
				  struct psc *psc, double *result)
{
  mparticles_t mprts = psc->particles->get_as<mparticles_t>();

  for (int p = 0; p < mprts.n_patches(); p++) {
    do_particle_energy(psc, mprts, p, result);
  }

  mprts.put_as(psc->particles, MP_DONT_COPY);
}

// ======================================================================
// psc_diag_item_particle_energy

struct psc_diag_item_ops psc_diag_item_particle_energy_ops = {
  .name      = "particle_energy",
  .run       = psc_diag_item_particle_energy_run,
  .nr_values = 2,
  .title     = { "E_electron", "E_ion" },
};

