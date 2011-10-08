
#include "psc.h"
#include "psc_particles_c.h"

#include <mrc_io.h>
#include <mrc_profile.h>
#include <stdlib.h>
#include <assert.h>

void
particles_c_alloc(particles_c_t *pp, int n_part)
{
  pp->n_alloced = n_part * 1.2;
  pp->particles = calloc(pp->n_alloced, sizeof(*pp->particles));
}

void
particles_c_realloc(particles_c_t *pp, int new_n_part)
{
  if (new_n_part <= pp->n_alloced)
    return;

  pp->n_alloced = new_n_part * 1.2;
  pp->particles = realloc(pp->particles, pp->n_alloced * sizeof(*pp->particles));
}

void
particles_c_free(particles_c_t *pp)
{
  free(pp->particles);
  pp->n_alloced = 0;
  pp->particles = NULL;
}

mparticles_c_t *
psc_mparticles_c_get_c(void *_particles_base)
{
  return _particles_base;
}

void
psc_mparticles_c_put_c(mparticles_c_t *particles, void *_particles_base)
{
}

static bool __gotten;

mparticles_c_t *
psc_mparticles_fortran_get_c(void *_particles_base)
{
  static int pr;
  if (!pr) {
    pr = prof_register("mparticles_get_c", 1., 0, 0);
  }
  prof_start(pr);

  assert(!__gotten);
  __gotten = true;
    
  mparticles_fortran_t *particles_base = _particles_base;

  int *nr_particles_by_patch = malloc(particles_base->nr_patches * sizeof(int));
  for (int p = 0; p < particles_base->nr_patches; p++) {
    nr_particles_by_patch[p] = psc_mparticles_fortran_nr_particles_by_patch(particles_base, p);
  }
  struct mrc_domain *domain = particles_base->domain;
  mparticles_c_t *particles = psc_mparticles_c_create(mrc_domain_comm(domain));
  psc_mparticles_c_set_domain_nr_particles(particles, domain, nr_particles_by_patch);
  psc_mparticles_c_setup(particles);
  free(nr_particles_by_patch);

  psc_foreach_patch(ppsc, p) {
    particles_fortran_t *pp_base = psc_mparticles_get_patch_fortran(particles_base, p);
    particles_c_t *pp = psc_mparticles_get_patch_c(particles, p);
    pp->n_part = pp_base->n_part;
    pp->particles = calloc(pp->n_part, sizeof(*pp->particles));
    for (int n = 0; n < pp_base->n_part; n++) {
      particle_fortran_t *part_base = particles_fortran_get_one(pp_base, n);
      particle_c_t *part = particles_c_get_one(pp, n);
      
      part->xi  = part_base->xi;
      part->yi  = part_base->yi;
      part->zi  = part_base->zi;
      part->pxi = part_base->pxi;
      part->pyi = part_base->pyi;
      part->pzi = part_base->pzi;
      part->qni = part_base->qni;
      part->mni = part_base->mni;
      part->wni = part_base->wni;
    }
  }

  prof_stop(pr);
  return particles;
}

void
psc_mparticles_fortran_put_c(mparticles_c_t *particles, void *_particles_base)
{
  static int pr;
  if (!pr) {
    pr = prof_register("mparticles_put_c", 1., 0, 0);
  }
  prof_start(pr);

  assert(__gotten);
  __gotten = false;

  mparticles_fortran_t *particles_base = _particles_base;
  psc_foreach_patch(ppsc, p) {
    particles_fortran_t *pp_base = psc_mparticles_get_patch_fortran(particles_base, p);
    particles_c_t *pp = psc_mparticles_get_patch_c(particles, p);
    assert(pp->n_part == pp_base->n_part);
    for (int n = 0; n < pp_base->n_part; n++) {
      particle_fortran_t *part_base = particles_fortran_get_one(pp_base, n);
      particle_c_t *part = particles_c_get_one(pp, n);
      
      part_base->xi  = part->xi;
      part_base->yi  = part->yi;
      part_base->zi  = part->zi;
      part_base->pxi = part->pxi;
      part_base->pyi = part->pyi;
      part_base->pzi = part->pzi;
      part_base->qni = part->qni;
      part_base->mni = part->mni;
      part_base->wni = part->wni;
    }
  }
  psc_mparticles_c_destroy(particles);

  prof_stop(pr);
}

// ======================================================================
// psc_mparticles_c

static void
_psc_mparticles_c_set_domain_nr_particles(mparticles_c_t *mparticles,
					  struct mrc_domain *domain,
					  int *nr_particles_by_patch)
{
  mparticles->domain = domain;
  mrc_domain_get_patches(domain, &mparticles->nr_patches);

  mparticles->data = calloc(mparticles->nr_patches, sizeof(*mparticles->data));
  for (int p = 0; p < mparticles->nr_patches; p++) {
    particles_c_alloc(psc_mparticles_get_patch_c(mparticles, p),
		      nr_particles_by_patch[p]);
  }
}

static void
_psc_mparticles_c_destroy(mparticles_c_t *mparticles)
{
  for (int p = 0; p < mparticles->nr_patches; p++) {
    particles_c_free(psc_mparticles_get_patch_c(mparticles, p));
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
  const char *path = psc_mparticles_c_name(mparticles);
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
  const char *path = psc_mparticles_c_name(mparticles);
  mrc_io_read_attr_int(io, path, "nr_patches", &mparticles->nr_patches);

  long h5_file;
  mrc_io_get_h5_file(io, &h5_file);
  hid_t group = H5Gopen(h5_file, path, H5P_DEFAULT); H5_CHK(group);
  mparticles->data = calloc(mparticles->nr_patches, sizeof(*mparticles->data));

  for (int p = 0; p < mparticles->nr_patches; p++) {
    particles_c_t *particles = psc_mparticles_get_patch_c(mparticles, p);
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

// ======================================================================
// psc_mparticles: subclass "c"
  
struct psc_mparticles_c_ops psc_mparticles_c_ops = {
  .name                    = "c",
  .set_domain_nr_particles = _psc_mparticles_c_set_domain_nr_particles,
};

static void
psc_mparticles_c_init()
{
  mrc_class_register_subclass(&mrc_class_psc_mparticles_c, &psc_mparticles_c_ops);
}

struct mrc_class_psc_mparticles_c mrc_class_psc_mparticles_c = {
  .name             = "psc_mparticles_c",
  .size             = sizeof(struct psc_mparticles_c),
  .init             = psc_mparticles_c_init,
  .destroy          = _psc_mparticles_c_destroy,
#ifdef HAVE_LIBHDF5_HL
  .write            = _psc_mparticles_c_write,
  .read             = _psc_mparticles_c_read,
#endif
};

