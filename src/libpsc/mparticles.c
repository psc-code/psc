
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
  mparticles->domain = domain;
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

#ifdef HAVE_LIBHDF5_HL

// FIXME. This is a rather bad break of proper layering, HDF5 should be all
// mrc_io business. OTOH, it could be called flexibility...

#include <hdf5.h>
#include <hdf5_hl.h>

#define H5_CHK(ierr) assert(ierr >= 0)
#define CE assert(ierr == 0)

static void
_psc_mparticles_c_write(mparticles_c_t *mparticles, struct mrc_io *io)
{
  int ierr;
  const char *path = psc_mparticles_c_name(mparticles);
  mrc_io_write_attr_int(io, path, "nr_patches", mparticles->nr_patches);
  
  long h5_file;
  mrc_io_get_h5_file(io, &h5_file);
  hid_t group = H5Gopen(h5_file, path, H5P_DEFAULT); H5_CHK(group);
  for (int p = 0; p < mparticles->nr_patches; p++) {
    particles_c_t *particles = &mparticles->p[p];
    char name[10]; sprintf(name, "p%d", p);

    hid_t groupp = H5Gcreate(group, name, H5P_DEFAULT, H5P_DEFAULT,
			     H5P_DEFAULT); H5_CHK(groupp);
    // save/restore n_alloced, too?
    ierr = H5LTset_attribute_int(groupp, ".", "n_part",
				 &particles->n_part, 1); CE;
    if (particles->n_part > 0) {
      hsize_t hdims[2] = { particles->n_part, 9 };
      ierr = H5LTmake_dataset_double(groupp, "particles_c", 2, hdims,
				     (double *) particles->particles); CE;
    }
    ierr = H5Gclose(groupp); CE;
  }

  ierr = H5Gclose(group); CE;
}

static void
_psc_mparticles_c_read(mparticles_c_t *mparticles, struct mrc_io *io)
{
  int ierr;
  const char *path = psc_mparticles_c_name(mparticles);
  mrc_io_read_attr_int(io, path, "nr_patches", &mparticles->nr_patches);

  long h5_file;
  mrc_io_get_h5_file(io, &h5_file);
  hid_t group = H5Gopen(h5_file, path, H5P_DEFAULT); H5_CHK(group);
  mparticles->p = calloc(mparticles->nr_patches, sizeof(*mparticles->p));

  for (int p = 0; p < mparticles->nr_patches; p++) {
    particles_c_t *particles = &mparticles->p[p];
    char name[10]; sprintf(name, "p%d", p);
    hid_t groupp = H5Gopen(group, name, H5P_DEFAULT); H5_CHK(groupp);
    int n_part;
    ierr = H5LTget_attribute_int(groupp, ".", "n_part", &n_part); CE;
    particles_c_alloc(particles, n_part);
    particles->n_part = n_part;
    if (n_part > 0) {
      ierr = H5LTread_dataset_double(groupp, "particles_c",
				     (double *) particles->particles); CE;
    }

    ierr = H5Gclose(groupp); CE;
  }

  ierr = H5Gclose(group); CE;
}

#endif


struct mrc_class_psc_mparticles_c mrc_class_psc_mparticles_c = {
  .name             = "psc_mparticles_c",
  .size             = sizeof(struct psc_mparticles_c),
  .destroy          = _psc_mparticles_c_destroy,
#ifdef HAVE_LIBHDF5_HL
  .write            = _psc_mparticles_c_write,
  .read             = _psc_mparticles_c_read,
#endif
};

