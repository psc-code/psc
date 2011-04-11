
#include "psc.h"

#include <stdlib.h>

// ----------------------------------------------------------------------
// mparticles_base_alloc

void
mparticles_base_alloc(struct mrc_domain *domain, mparticles_base_t *particles,
		      int *nr_particles_by_patch)
{
  mrc_domain_get_patches(domain, &particles->nr_patches);

  particles->p = calloc(particles->nr_patches, sizeof(*particles->p));
  for (int p = 0; p < particles->nr_patches; p++) {
    particles_base_alloc(&particles->p[p], nr_particles_by_patch[p]);
  }
}

// ----------------------------------------------------------------------
// mparticles_base_destroy

void
mparticles_base_destroy(mparticles_base_t *particles)
{
  for (int p = 0; p < particles->nr_patches; p++) {
    particles_base_free(&particles->p[p]);
  }
  free(particles->p);
  particles->nr_patches = -1;
}

