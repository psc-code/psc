
#include "psc.h"
#include "psc_particles_c.h"

#include <mrc_io.h>
#include <stdlib.h>
#include <assert.h>

static void
_psc_mparticles_c_alloc_patch(mparticles_c_t *mp, int p, int n_part)
{
  particles_c_t *pp = psc_mparticles_get_patch_c(mp, p);
  pp->n_part = n_part;
  pp->n_alloced = n_part * 1.2;
  pp->particles = calloc(pp->n_alloced, sizeof(*pp->particles));
}

static void
_psc_mparticles_c_free_patch(mparticles_c_t *mp, int p)
{
  particles_c_t *pp = psc_mparticles_get_patch_c(mp, p);
  free(pp->particles);
  pp->n_alloced = 0;
  pp->particles = NULL;
}

void
particles_c_realloc(particles_c_t *pp, int new_n_part)
{
  if (new_n_part <= pp->n_alloced)
    return;

  pp->n_alloced = new_n_part * 1.2;
  pp->particles = realloc(pp->particles, pp->n_alloced * sizeof(*pp->particles));
}

// ======================================================================
// psc_mparticles_c

static void
_psc_mparticles_c_setup(mparticles_c_t *mparticles)
{
  assert(mparticles->nr_particles_by_patch);

  mparticles->data = calloc(mparticles->nr_patches, sizeof(particles_c_t));
  for (int p = 0; p < mparticles->nr_patches; p++) {
    _psc_mparticles_c_alloc_patch(mparticles, p, mparticles->nr_particles_by_patch[p]);
  }

  free(mparticles->nr_particles_by_patch);
  mparticles->nr_particles_by_patch = NULL;
}
									
static int
_psc_mparticles_c_nr_particles_by_patch(mparticles_c_t *mparticles, int p)
{
  return psc_mparticles_get_patch_c(mparticles, p)->n_part;
}

static void
_psc_mparticles_c_destroy(mparticles_c_t *mparticles)
{
  for (int p = 0; p < mparticles->nr_patches; p++) {
    _psc_mparticles_c_free_patch(mparticles, p);
  }
  free(mparticles->data);
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
  const char *path = psc_mparticles_name(mparticles);
  mrc_io_write_attr_int(io, path, "nr_patches", mparticles->nr_patches);
  
  long h5_file;
  mrc_io_get_h5_file(io, &h5_file);
  hid_t group = H5Gopen(h5_file, path, H5P_DEFAULT); H5_CHK(group);
  for (int p = 0; p < mparticles->nr_patches; p++) {
    particles_c_t *particles = psc_mparticles_get_patch_c(mparticles, p);
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
  const char *path = psc_mparticles_name(mparticles);
  mrc_io_read_attr_int(io, path, "nr_patches", &mparticles->nr_patches);

  long h5_file;
  mrc_io_get_h5_file(io, &h5_file);
  hid_t group = H5Gopen(h5_file, path, H5P_DEFAULT); H5_CHK(group);
  mparticles->data = calloc(mparticles->nr_patches, sizeof(particles_c_t));

  for (int p = 0; p < mparticles->nr_patches; p++) {
    particles_c_t *particles = psc_mparticles_get_patch_c(mparticles, p);
    char name[10]; sprintf(name, "p%d", p);
    hid_t groupp = H5Gopen(group, name, H5P_DEFAULT); H5_CHK(groupp);
    int n_part;
    ierr = H5LTget_attribute_int(groupp, ".", "n_part", &n_part); CE;
    _psc_mparticles_c_alloc_patch(mparticles, p, n_part);
    if (n_part > 0) {
      ierr = H5LTread_dataset_double(groupp, "particles_c",
				     (double *) particles->particles); CE;
    }

    ierr = H5Gclose(groupp); CE;
  }

  ierr = H5Gclose(group); CE;
}

#endif

// ======================================================================
// psc_mparticles: subclass "c"
  
struct psc_mparticles_ops psc_mparticles_c_ops = {
  .name                    = "c",
  .setup                   = _psc_mparticles_c_setup,
  .destroy                 = _psc_mparticles_c_destroy,
#ifdef HAVE_LIBHDF5_HL
  .write                   = _psc_mparticles_c_write,
  .read                    = _psc_mparticles_c_read,
#endif
  .nr_particles_by_patch   = _psc_mparticles_c_nr_particles_by_patch,
};

