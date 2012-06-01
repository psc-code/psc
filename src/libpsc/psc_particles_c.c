
#include "psc.h"
#include "psc_particles_c.h"

#include <mrc_io.h>
#include <stdlib.h>
#include <assert.h>

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

#ifdef HAVE_LIBHDF5_HL

// FIXME. This is a rather bad break of proper layering, HDF5 should be all
// mrc_io business. OTOH, it could be called flexibility...

#include <hdf5.h>
#include <hdf5_hl.h>

#define H5_CHK(ierr) assert(ierr >= 0)
#define CE assert(ierr == 0)

// ----------------------------------------------------------------------
// psc_particles_c_write

static void
psc_particles_c_write(struct psc_particles *prts, struct mrc_io *io)
{
  int ierr;
  assert(sizeof(particle_c_t) / sizeof(particle_c_real_t) == 10);
  assert(sizeof(particle_c_real_t) == sizeof(double));

  long h5_file;
  mrc_io_get_h5_file(io, &h5_file);

  hid_t group = H5Gopen(h5_file, psc_particles_name(prts), H5P_DEFAULT); H5_CHK(group);
  // save/restore n_alloced, too?
  ierr = H5LTset_attribute_int(group, ".", "p", &prts->p, 1); CE;
  ierr = H5LTset_attribute_int(group, ".", "n_part", &prts->n_part, 1); CE;
  ierr = H5LTset_attribute_uint(group, ".", "flags", &prts->flags, 1); CE;
  if (prts->n_part > 0) {
    // in a rather ugly way, we write the long "kind" member as a double
    hsize_t hdims[2] = { prts->n_part, 10 };
    ierr = H5LTmake_dataset_double(group, "particles_c", 2, hdims,
				   (double *) particles_c_get_one(prts, 0)); CE;
  }
  ierr = H5Gclose(group); CE;
}

// ----------------------------------------------------------------------
// psc_particles_c_read

static void
psc_particles_c_read(struct psc_particles *prts, struct mrc_io *io)
{
  int ierr;
  long h5_file;
  mrc_io_get_h5_file(io, &h5_file);

  hid_t group = H5Gopen(h5_file, psc_particles_name(prts), H5P_DEFAULT); H5_CHK(group);
  ierr = H5LTget_attribute_int(group, ".", "p", &prts->p); CE;
  ierr = H5LTget_attribute_int(group, ".", "n_part", &prts->n_part); CE;
  ierr = H5LTget_attribute_uint(group, ".", "flags", &prts->flags); CE;
  psc_particles_setup(prts);
  if (prts->n_part > 0) {
    ierr = H5LTread_dataset_double(group, "particles_c",
				   (double *) particles_c_get_one(prts, 0)); CE;
  }
  ierr = H5Gclose(group); CE;
}

#endif

// ======================================================================
// psc_mparticles_c

static void
_psc_mparticles_c_write(struct psc_mparticles *mparticles, struct mrc_io *io)
{
  const char *path = psc_mparticles_name(mparticles);
  mrc_io_write_obj_ref(io, path, "domain", (struct mrc_obj *) mparticles->domain);
  mrc_io_write_attr_int(io, path, "nr_patches", mparticles->nr_patches);
  mrc_io_write_attr_int(io, path, "flags", mparticles->flags);
  
  for (int p = 0; p < mparticles->nr_patches; p++) {
    psc_particles_write(mparticles->prts[p], io);
  }
}

static void
_psc_mparticles_c_read(mparticles_c_t *mparticles, struct mrc_io *io)
{
  const char *path = psc_mparticles_name(mparticles);
  mparticles->domain = (struct mrc_domain *)
    mrc_io_read_obj_ref(io, path, "domain", &mrc_class_mrc_domain);
  mrc_io_read_attr_int(io, path, "nr_patches", &mparticles->nr_patches);
  mrc_io_read_attr_int(io, path, "flags", (int *) &mparticles->flags);

  mparticles->prts = calloc(mparticles->nr_patches, sizeof(*mparticles->prts));
  mparticles->nr_particles_by_patch =
    calloc(mparticles->nr_patches, sizeof(*mparticles->nr_particles_by_patch));
  for (int p = 0; p < mparticles->nr_patches; p++) {
    char name[20]; sprintf(name, "prts%d", p);
    mparticles->prts[p] = psc_particles_read(io, name);
  }
}

// ======================================================================
// psc_mparticles: subclass "c"
  
struct psc_mparticles_ops psc_mparticles_c_ops = {
  .name                    = "c",
  .write                   = _psc_mparticles_c_write,
  .read                    = _psc_mparticles_c_read,
};

// ======================================================================
// psc_particles: subclass "c"

struct psc_particles_ops psc_particles_c_ops = {
  .name                    = "c",
  .size                    = sizeof(struct psc_particles_c),
  .setup                   = psc_particles_c_setup,
  .destroy                 = psc_particles_c_destroy,
#ifdef HAVE_LIBHDF5_HL
  .write                   = psc_particles_c_write,
  .read                    = psc_particles_c_read,
#endif
};
