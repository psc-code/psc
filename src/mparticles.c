
#include "psc.h"

// ----------------------------------------------------------------------
// mparticles_base_destroy

void
mparticles_base_destroy(mparticles_base_t *particles)
{
  foreach_patch(p) {
    particles_base_free(&particles->p[p]);
  }
  free(particles->p);
}

