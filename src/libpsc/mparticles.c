
#include "psc.h"

#include <stdlib.h>

#define MAKE_MPARTICLES_METHODS(type)					\
									\
mparticles_##type##_t *							\
mparticles_##type##_create(MPI_Comm comm)				\
{									\
  mparticles_##type##_t *mparticles = calloc(1, sizeof(*mparticles));	\
									\
  return mparticles;							\
}									\
									\
void									\
mparticles_##type##_alloc(mparticles_##type##_t *mparticles,		\
			  struct mrc_domain *domain,			\
			  int *nr_particles_by_patch)			\
{									\
  mrc_domain_get_patches(domain, &mparticles->nr_patches);		\
									\
  mparticles->p = calloc(mparticles->nr_patches,			\
			 sizeof(*mparticles->p));			\
  for (int p = 0; p < mparticles->nr_patches; p++) {			\
    particles_##type##_alloc(&mparticles->p[p],				\
			     nr_particles_by_patch[p]);			\
  }									\
									\
  return mparticles;							\
}									\
									\
void									\
mparticles_##type##_destroy(mparticles_##type##_t *particles)		\
{									\
  for (int p = 0; p < particles->nr_patches; p++) {			\
    particles_##type##_free(&particles->p[p]);				\
  }									\
  free(particles->p);							\
  free(particles);							\
}

#if PARTICLES_BASE == PARTICLES_FORTRAN
MAKE_MPARTICLES_METHODS(fortran)
#elif PARTICLES_BASE == PARTICLES_C
MAKE_MPARTICLES_METHODS(c)
#elif PARTICLES_BASE == PARTICLES_SSE2
MAKE_MPARTICLES_METHODS(sse2)
#elif PARTICLES_BASE == PARTICLES_CBE
MAKE_MPARTICLES_METHODS(cbe)
#endif
