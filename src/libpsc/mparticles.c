
#include "psc.h"

#include <stdlib.h>

#define MAKE_MPARTICLES_METHODS(type)					\
									\
void									\
psc_mparticles_##type##_set_domain_nr_particles(mparticles_##type##_t *mparticles,		\
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
static void								\
_psc_mparticles_##type##_destroy(mparticles_##type##_t *mparticles)	\
{									\
  for (int p = 0; p < mparticles->nr_patches; p++) {			\
    particles_##type##_free(&mparticles->p[p]);				\
  }									\
  free(mparticles->p);							\
}									\
									\
struct mrc_class_psc_mparticles_##type mrc_class_psc_mparticles_##type = {	\
  .name             = "psc_mparticles_" #type,				\
  .size             = sizeof(struct psc_mparticles_##type),		\
  .destroy          = _psc_mparticles_##type##_destroy,			\
};


#if PARTICLES_BASE == PARTICLES_FORTRAN
MAKE_MPARTICLES_METHODS(fortran)
#elif PARTICLES_BASE == PARTICLES_SSE2
MAKE_MPARTICLES_METHODS(sse2)
#elif PARTICLES_BASE == PARTICLES_CBE
MAKE_MPARTICLES_METHODS(cbe)
#endif

// ======================================================================
// psc_mparticles_c

#include <mrc_io.h>

void
psc_mparticles_c_set_domain_nr_particles(mparticles_c_t *mparticles,
					    struct mrc_domain *domain,
					    int *nr_particles_by_patch)
{
  mrc_domain_get_patches(domain, &mparticles->nr_patches);

  mparticles->p = calloc(mparticles->nr_patches, sizeof(*mparticles->p));
  for (int p = 0; p < mparticles->nr_patches; p++) {
    particles_c_alloc(&mparticles->p[p],
		      nr_particles_by_patch[p]);
  }
}

static void
_psc_mparticles_c_destroy(mparticles_c_t *mparticles)
{
  for (int p = 0; p < mparticles->nr_patches; p++) {
    particles_c_free(&mparticles->p[p]);
  }
  free(mparticles->p);
}

static void
_psc_mparticles_c_write(mparticles_c_t *mparticles, struct mrc_io *io)
{
  const char *path = psc_mparticles_c_name(mparticles);
  mrc_io_write_attr_int(io, path, "nr_patches", mparticles->nr_patches);
  
  for (int p = 0; p < mparticles->nr_patches; p++) {
    char name[10]; sprintf(name, "p%d", p);
    // FIXME, should use n_part, maybe n_alloced, too
    mrc_io_write_attr_int(io, path, name, mparticles->p[p].n_alloced);
  }
}

static void
_psc_mparticles_c_read(mparticles_c_t *mparticles, struct mrc_io *io)
{
  const char *path = psc_mparticles_c_name(mparticles);
  mrc_io_read_attr_int(io, path, "nr_patches", &mparticles->nr_patches);

  int nr_particles_by_patch[mparticles->nr_patches];
  for (int p = 0; p < mparticles->nr_patches; p++) {
    char name[10]; sprintf(name, "p%d", p);
    mrc_io_read_attr_int(io, path, name, &nr_particles_by_patch[p]);
  }

  // this should probably move into ->setup(), joined with 
  // the 2nd half of set_domain_nr_particles...()
  mparticles->p = calloc(mparticles->nr_patches, sizeof(*mparticles->p));
  for (int p = 0; p < mparticles->nr_patches; p++) {
    particles_c_alloc(&mparticles->p[p], nr_particles_by_patch[p]);
  }
}

struct mrc_class_psc_mparticles_c mrc_class_psc_mparticles_c = {
  .name             = "psc_mparticles_c",
  .size             = sizeof(struct psc_mparticles_c),
  .destroy          = _psc_mparticles_c_destroy,
  .write            = _psc_mparticles_c_write,
  .read             = _psc_mparticles_c_read,
};

