
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
mparticles_##type##_set_domain_nr_particles(mparticles_##type##_t *mparticles,		\
					    struct mrc_domain *domain,	\
					    int *nr_particles_by_patch)	\
{									\
  mrc_domain_get_patches(domain, &mparticles->nr_patches);		\
									\
  mparticles->p = calloc(mparticles->nr_patches,			\
			 sizeof(*mparticles->p));			\
  for (int p = 0; p < mparticles->nr_patches; p++) {			\
    particles_##type##_alloc(&mparticles->p[p],				\
			     nr_particles_by_patch[p]);			\
  }									\
}									\
									\
void									\
mparticles_##type##_setup(mparticles_##type##_t *mparticles)		\
{									\
}									\
									\
void									\
mparticles_##type##_destroy(mparticles_##type##_t *mparticles)		\
{									\
  for (int p = 0; p < mparticles->nr_patches; p++) {			\
    particles_##type##_free(&mparticles->p[p]);				\
  }									\
  free(mparticles->p);							\
  free(mparticles);							\
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
