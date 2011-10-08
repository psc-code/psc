
#include "psc.h"

#include <stdlib.h>

#define MAKE_MPARTICLES_METHODS(type)					\
									\
static void								\
_psc_mparticles_##type##_set_domain_nr_particles(mparticles_##type##_t *mparticles, \
						 struct mrc_domain *domain, \
						 int *nr_particles_by_patch) \
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
struct psc_mparticles_##type##_ops psc_mparticles_##type##_ops = {	\
  .name                    = #type,					\
  .set_domain_nr_particles = _psc_mparticles_##type##_set_domain_nr_particles, \
};									\
									\
static void								\
psc_mparticles_##type##_init()						\
{									\
  mrc_class_register_subclass(&mrc_class_psc_mparticles_##type, &psc_mparticles_##type##_ops); \
}									\
									\
struct mrc_class_psc_mparticles_##type mrc_class_psc_mparticles_##type = {	\
  .name             = "psc_mparticles_" #type,				\
  .size             = sizeof(struct psc_mparticles_##type),		\
  .init		    = psc_mparticles_##type##_init,			\
  .destroy          = _psc_mparticles_##type##_destroy,			\
};									\


#define _MAKE_MPARTICLES_METHODS(type)					\
void									\
psc_mparticles_##type##_set_domain_nr_particles(mparticles_##type##_t *mparticles,	\
						struct mrc_domain *domain, \
						int *nr_particles_by_patch) \
{									\
  psc_mparticles_set_domain_nr_particles((struct psc_mparticles *) mparticles, \
					 domain, nr_particles_by_patch); \
}

MAKE_MPARTICLES_METHODS(fortran)
_MAKE_MPARTICLES_METHODS(fortran)
#ifdef USE_SSE2
MAKE_MPARTICLES_METHODS(sse2)
_MAKE_MPARTICLES_METHODS(sse2)
#endif
#ifdef USE_CBE
MAKE_MPARTICLES_METHODS(cbe)
_MAKE_MPARTICLES_METHODS(cbe)
#endif

_MAKE_MPARTICLES_METHODS(c)
_MAKE_MPARTICLES_METHODS(cuda)

// ======================================================================

void
psc_mparticles_set_domain_nr_particles(struct psc_mparticles *mparticles,
				       struct mrc_domain *domain,
				       int *nr_particles_by_patch)
{
  struct psc_mparticles_ops *ops = psc_mparticles_ops(mparticles);
  assert(ops && ops->set_domain_nr_particles);
  return ops->set_domain_nr_particles(mparticles, domain, nr_particles_by_patch);
}


