
#include "psc.h"

#include <stdlib.h>

// ----------------------------------------------------------------------
// mparticles_base_alloc

void
mparticles_base_alloc(mparticles_base_t *particles, int *nr_particles_by_patch)
{
  particles->p = calloc(psc.nr_patches, sizeof(*particles->p));
  foreach_patch(p) {
    particles_base_alloc(&particles->p[p], nr_particles_by_patch[p]);
  }
}

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

