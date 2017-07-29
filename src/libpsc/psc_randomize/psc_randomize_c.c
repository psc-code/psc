
#include "psc_randomize_private.h"
#include "psc_particles_as_c.h"

#include <mrc_profile.h>

// FIXME, this is very inefficient if followed by a sort..

// ----------------------------------------------------------------------
// psc_randomize_c_run

static void
psc_randomize_c_run(struct psc_randomize *randomize,
		    struct psc_mparticles *mprts_base)
{
  struct psc_mparticles *mprts = psc_mparticles_get_as(mprts_base, PARTICLE_TYPE, 0);

  for (int p = 0; p < mprts->nr_patches; p++) {
    struct psc_particles *prts = psc_mparticles_get_patch(mprts, p);
    for (int i = 0; i < prts->n_part; i++) {
      int j = random() % prts->n_part;
      particle_t tmp = *particles_get_one(prts, i);
      *particles_get_one(prts, i) = *particles_get_one(prts, j);
      *particles_get_one(prts, j) = tmp;
    }
  }

  psc_mparticles_put_as(mprts, mprts_base, 0);
}

// ======================================================================
// psc_randomize: subclass "c"

struct psc_randomize_ops psc_randomize_c_ops = {
  .name                  = "c",
  .run                   = psc_randomize_c_run,
};
