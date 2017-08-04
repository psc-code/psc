
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

  prts->n_alloced = psc_particles_size(prts) * 1.2;
  c->particles = calloc(prts->n_alloced, sizeof(*c->particles));
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
  if (new_n_part <= prts->n_alloced)
    return;

  prts->n_alloced = new_n_part * 1.2;
  c->particles = realloc(c->particles, prts->n_alloced * sizeof(*c->particles));
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

  hid_t group = H5Gopen(h5_file, mrc_io_obj_path(io, prts), H5P_DEFAULT); H5_CHK(group);
  // save/restore n_alloced, too?
  ierr = H5LTset_attribute_int(group, ".", "p", &prts->p, 1); CE;
  int n_prts = psc_particles_size(prts);
  ierr = H5LTset_attribute_int(group, ".", "n_part", &n_prts, 1); CE;
  ierr = H5LTset_attribute_uint(group, ".", "flags", &prts->flags, 1); CE;
  if (n_prts > 0) {
    // in a rather ugly way, we write the long "kind" member as a double
    hsize_t hdims[2] = { n_prts, 10 };
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

  hid_t group = H5Gopen(h5_file, mrc_io_obj_path(io, prts), H5P_DEFAULT); H5_CHK(group);
  ierr = H5LTget_attribute_int(group, ".", "p", &prts->p); CE;
  int n_prts;
  ierr = H5LTget_attribute_int(group, ".", "n_part", &n_prts); CE;
  psc_particles_resize(prts, n_prts);
  ierr = H5LTget_attribute_uint(group, ".", "flags", &prts->flags); CE;
  psc_particles_setup(prts);
  if (n_prts > 0) {
    ierr = H5LTread_dataset_double(group, "particles_c",
				   (double *) particles_c_get_one(prts, 0)); CE;
  }
  ierr = H5Gclose(group); CE;
}

#endif

// ======================================================================
// psc_mparticles: subclass "c"
  
struct psc_mparticles_ops psc_mparticles_c_ops = {
  .name                    = "c",
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
