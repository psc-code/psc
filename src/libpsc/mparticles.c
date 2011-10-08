
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
static int									\
_psc_mparticles_##type##_nr_particles_by_patch(mparticles_##type##_t *mparticles, \
					      int p)			\
{									\
  return psc_mparticles_get_patch_##type(mparticles, p)->n_part;	\
}									\
									\
struct psc_mparticles_##type##_ops psc_mparticles_##type##_ops = {	\
  .name                    = #type,					\
  .set_domain_nr_particles = _psc_mparticles_##type##_set_domain_nr_particles, \
  .nr_particles_by_patch = _psc_mparticles_##type##_nr_particles_by_patch, \
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
}									\
									\
int									\
psc_mparticles_##type##_nr_particles_by_patch(mparticles_##type##_t *mparticles, \
					      int p) \
{									\
  return psc_mparticles_nr_particles_by_patch((struct psc_mparticles *) mparticles, p); \
}									\
									\
mparticles_c_t *							\
psc_mparticles_##type##_get_c(void *particles_base)			\
{									\
  return psc_mparticles_get_c(particles_base);				\
}									\
									\
void									\
psc_mparticles_##type##_put_c(mparticles_c_t *particles, void *particles_base) \
{									\
  psc_mparticles_put_c(particles, particles_base);			\
}									\


#ifdef USE_SSE2
MAKE_MPARTICLES_METHODS(sse2)
_MAKE_MPARTICLES_METHODS(sse2)
#endif
#ifdef USE_CBE
MAKE_MPARTICLES_METHODS(cbe)
_MAKE_MPARTICLES_METHODS(cbe)
#endif

_MAKE_MPARTICLES_METHODS(fortran)
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

int
psc_mparticles_nr_particles_by_patch(struct psc_mparticles *mparticles, int p)
{
  struct psc_mparticles_ops *ops = psc_mparticles_ops(mparticles);
  assert(ops && ops->nr_particles_by_patch);
  return ops->nr_particles_by_patch(mparticles, p);
}

mparticles_c_t *
psc_mparticles_get_c(struct psc_mparticles *particles_base)
{
  struct psc_mparticles_ops *ops = psc_mparticles_ops(particles_base);
  assert(ops && ops->get_c);
  return (mparticles_c_t *) ops->get_c(particles_base);
}

void
psc_mparticles_put_c(mparticles_c_t *particles, struct psc_mparticles *particles_base)
{
  struct psc_mparticles_ops *ops = psc_mparticles_ops(particles_base);
  assert(ops && ops->put_c);
  ops->put_c((struct psc_mparticles *) particles, particles_base);
}



