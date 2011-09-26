
#include "psc_randomize_private.h"

#include <mrc_profile.h>

// FIXME, this is very inefficient if followed by a sort..

// ----------------------------------------------------------------------
// psc_randomize_c_run

static void
psc_randomize_c_run(struct psc_randomize *randomize,
		    mparticles_base_t *particles_base)
{
  psc_foreach_patch(ppsc, p) {
    particles_base_t *pp = &particles_base->p[p];
    for (int i = 0; i < pp->n_part; i++) {
      int j = random() % pp->n_part;
      particle_base_t tmp = pp->particles[i];
      pp->particles[i] = pp->particles[j];
      pp->particles[j] = tmp;
    }
  }
}

// ======================================================================
// psc_randomize: subclass "c"

struct psc_randomize_ops psc_randomize_c_ops = {
  .name                  = "c",
  .run                   = psc_randomize_c_run,
};
