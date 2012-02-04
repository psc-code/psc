
#include "psc.h"

#include <mrc_profile.h>
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
			    sizeof(particles_##type##_t));		\
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


#ifdef USE_SSE2
MAKE_MPARTICLES_METHODS(sse2)
#endif
#ifdef USE_CBE
MAKE_MPARTICLES_METHODS(cbe)
#endif

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

#define MAKE_MPARTICLES_GET_PUT(type)					\
									\
mparticles_##type##_t *							\
psc_mparticles_get_##type(struct psc_mparticles *particles_base,	\
			  unsigned int flags)				\
{									\
  struct psc_mparticles_ops *ops = psc_mparticles_ops(particles_base);	\
  if (ops == &psc_mparticles_##type##_ops) {				\
    return particles_base;						\
  }									\
  static int pr;							\
  if (!pr) {								\
    pr = prof_register("mparticles_get_" #type, 1., 0, 0);		\
  }									\
  prof_start(pr);							\
  assert(ops);								\
  if (!ops->copy_to_##type) {						\
    fprintf(stderr, "ERROR: missing copy_to_" #type			\
	    " in psc_mparticles '%s'!\n",				\
	    psc_mparticles_type(particles_base));			\
    assert(0);								\
  }									\
									\
  int *nr_particles_by_patch =						\
    malloc(particles_base->nr_patches * sizeof(int));			\
  for (int p = 0; p < particles_base->nr_patches; p++) {		\
  nr_particles_by_patch[p] =						\
    psc_mparticles_nr_particles_by_patch(particles_base, p);		\
  }									\
  struct psc_mparticles *mp =						\
    psc_mparticles_create(psc_mparticles_comm(particles_base));		\
  psc_mparticles_set_type(mp, #type);					\
  psc_mparticles_set_domain_nr_particles(mp,				\
					 particles_base->domain,	\
					 nr_particles_by_patch);	\
  psc_mparticles_setup(mp);						\
  free(nr_particles_by_patch);						\
									\
  ops->copy_to_##type(particles_base, mp, flags);			\
  prof_stop(pr);							\
  return mp;								\
}									\
									\
void									\
psc_mparticles_put_##type(mparticles_##type##_t *particles,		\
		     struct psc_mparticles *particles_base)		\
{									\
  struct psc_mparticles_ops *ops = psc_mparticles_ops(particles_base);	\
  if (ops == &psc_mparticles_##type##_ops) {				\
    return;								\
  }									\
  static int pr;							\
  if (!pr) {								\
    pr = prof_register("mparticles_get_" #type, 1., 0, 0);		\
  }									\
  prof_start(pr);							\
  assert(ops);								\
  if (!ops->copy_from_##type) {						\
    fprintf(stderr, "ERROR: missing copy_from_"#type			\
	    " in psc_mparticles '%s'!\n",				\
	    psc_mparticles_type(particles_base));			\
    assert(0);								\
  }									\
  ops->copy_from_##type(particles_base, particles,			\
			MP_NEED_BLOCK_OFFSETS | MP_NEED_CELL_OFFSETS);	\
  psc_mparticles_destroy(particles);					\
  prof_stop(pr);							\
}									\

MAKE_MPARTICLES_GET_PUT(c)
MAKE_MPARTICLES_GET_PUT(fortran)
#ifdef USE_CUDA
MAKE_MPARTICLES_GET_PUT(cuda)
#endif

// ======================================================================

static void
psc_mparticles_init()
{
  mrc_class_register_subclass(&mrc_class_psc_mparticles, &psc_mparticles_c_ops);
  mrc_class_register_subclass(&mrc_class_psc_mparticles, &psc_mparticles_fortran_ops);
#ifdef USE_CUDA
  mrc_class_register_subclass(&mrc_class_psc_mparticles, &psc_mparticles_cuda_ops);
#endif
}

struct mrc_class_psc_mparticles mrc_class_psc_mparticles = {
  .name             = "psc_mparticles_c",
  .size             = sizeof(struct psc_mparticles),
  .init             = psc_mparticles_init,
};

