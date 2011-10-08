
#include "psc.h"

#include <stdlib.h>

#define MAKE_MPARTICLES_METHODS(type)					\
									\
void									\
psc_mparticles_##type##_set_domain_nr_particles(mparticles_##type##_t *mparticles,		\
					    struct mrc_domain *domain,	\
					    int *nr_particles_by_patch)	\
{									\
  mparticles->domain = domain;						\
  mrc_domain_get_patches(domain, &mparticles->nr_patches);		\
									\
  mparticles->data = calloc(mparticles->nr_patches,			\
			    sizeof(*mparticles->data));			\
  for (int p = 0; p < mparticles->nr_patches; p++) {			\
    particles_##type##_alloc(&mparticles->data[p],			\
			     nr_particles_by_patch[p]);			\
  }									\
}									\
									\
static void								\
_psc_mparticles_##type##_destroy(mparticles_##type##_t *mparticles)	\
{									\
  for (int p = 0; p < mparticles->nr_patches; p++) {			\
    particles_##type##_free(&mparticles->data[p]);			\
  }									\
  free(mparticles->data);						\
}									\
									\
struct mrc_class_psc_mparticles_##type mrc_class_psc_mparticles_##type = {	\
  .name             = "psc_mparticles_" #type,				\
  .size             = sizeof(struct psc_mparticles_##type),		\
  .destroy          = _psc_mparticles_##type##_destroy,			\
};


MAKE_MPARTICLES_METHODS(fortran)
#ifdef USE_SSE2
MAKE_MPARTICLES_METHODS(sse2)
#endif
#ifdef USE_CBE
MAKE_MPARTICLES_METHODS(cbe)
#endif

