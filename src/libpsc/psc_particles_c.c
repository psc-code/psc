
#include "psc.h"
#include "psc_particles_c.h"

#include <mrc_io.h>
#include <stdlib.h>
#include <assert.h>

static void *
_psc_mparticles_c_alloc_patch(int p, int n_part, unsigned int flags)
{
  MPI_Comm comm = MPI_COMM_WORLD; // FIXME!
  struct psc_particles *prts = psc_particles_create(comm);
  psc_particles_set_type(prts, "c");
  prts->n_part = n_part;
  psc_particles_setup(prts);
  return prts;
}

static void
_psc_mparticles_c_free_patch(int p, void *_pp)
{
  struct psc_particles *prts = _pp;
  psc_particles_destroy(prts);
}

// ======================================================================
// psc_particles "c"

static void
psc_particles_c_setup(struct psc_particles *prts)
{
  struct psc_particles_c *c = psc_particles_c(prts);

  c->n_alloced = prts->n_part * 1.2;
  c->particles = calloc(c->n_alloced, sizeof(*c->particles));
}

static void
psc_particles_c_destroy(struct psc_particles *prts)
{
  struct psc_particles_c *c = psc_particles_c(prts);

  free(c->particles);
}

void
particles_c_realloc(struct psc_particles *prts, int new_n_part)
{
  struct psc_particles_c *c = psc_particles_c(prts);
  if (new_n_part <= c->n_alloced)
    return;

  c->n_alloced = new_n_part * 1.2;
  c->particles = realloc(c->particles, c->n_alloced * sizeof(*c->particles));
}

// ======================================================================
// psc_mparticles_c

static int
_psc_mparticles_c_nr_particles_by_patch(mparticles_c_t *mparticles, int p)
{
  return psc_mparticles_get_patch(mparticles, p)->n_part;
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
    struct psc_particles *prts = psc_mparticles_get_patch(mparticles, p);
    struct psc_particles_c *c = psc_particles_c(prts);
    char name[10]; sprintf(name, "p%d", p);

    hid_t groupp = H5Gcreate(group, name, H5P_DEFAULT, H5P_DEFAULT,
			     H5P_DEFAULT); H5_CHK(groupp);
    // save/restore n_alloced, too?
    ierr = H5LTset_attribute_int(groupp, ".", "n_part",
				 &prts->n_part, 1); CE;
    if (prts->n_part > 0) {
      assert(0); // need to fix for "kind"
      hsize_t hdims[2] = { prts->n_part, 9 };
      ierr = H5LTmake_dataset_double(groupp, "particles_c", 2, hdims,
				     (double *) c->particles); CE;
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
  mparticles->patches = calloc(mparticles->nr_patches, sizeof(*mparticles->patches));
  mparticles->nr_particles_by_patch =
    calloc(mparticles->nr_patches, sizeof(*mparticles->nr_particles_by_patch));
  for (int p = 0; p < mparticles->nr_patches; p++) {
    char name[10]; sprintf(name, "p%d", p);
    hid_t groupp = H5Gopen(group, name, H5P_DEFAULT); H5_CHK(groupp);
    int n_part;
    ierr = H5LTget_attribute_int(groupp, ".", "n_part", &n_part); CE;
    struct psc_particles *prts = _psc_mparticles_c_alloc_patch(p, n_part, 0);
    struct psc_particles_c *c = psc_particles_c(prts);
    mparticles->patches[p] = prts;
    if (n_part > 0) {
      ierr = H5LTread_dataset_double(groupp, "particles_c",
				     (double *) c->particles); CE;
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
#ifdef HAVE_LIBHDF5_HL
  .write                   = _psc_mparticles_c_write,
  .read                    = _psc_mparticles_c_read,
#endif
  .nr_particles_by_patch   = _psc_mparticles_c_nr_particles_by_patch,
  .alloc_patch             = _psc_mparticles_c_alloc_patch,
  .free_patch              = _psc_mparticles_c_free_patch,
};

// ======================================================================
// psc_particles: subclass "c"

struct psc_particles_ops psc_particles_c_ops = {
  .name                    = "c",
  .size                    = sizeof(struct psc_particles_c),
  .setup                   = psc_particles_c_setup,
  .destroy                 = psc_particles_c_destroy,
};
