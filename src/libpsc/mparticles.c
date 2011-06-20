
#include "psc.h"

#include <stdlib.h>

// ----------------------------------------------------------------------
// mparticles_base_alloc

mparticles_base_t *
mparticles_base_alloc(struct mrc_domain *domain, int *nr_particles_by_patch)
{
  mparticles_base_t *mparticles = calloc(1, sizeof(*mparticles));

  mrc_domain_get_patches(domain, &mparticles->nr_patches);

  mparticles->p = calloc(mparticles->nr_patches, sizeof(*mparticles->p));
  for (int p = 0; p < mparticles->nr_patches; p++) {
    particles_base_alloc(&mparticles->p[p], nr_particles_by_patch[p]);
  }
  
  return mparticles;
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
  free(particles);
}

