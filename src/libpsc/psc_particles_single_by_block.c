
#include "psc_particles_single_by_block.h"

#include "psc.h"
#include "psc_particles_single.h"

#include <mrc_io.h>
#include <stdlib.h>
#include <assert.h>
#include <math.h>

// ======================================================================
// psc_particles "single_by_block"

static void
psc_particles_single_by_block_setup(struct psc_particles *prts)
{
  struct psc_particles_single_by_block *sub = psc_particles_single_by_block(prts);

  sub->n_alloced = prts->n_part * 1.2;
  sub->particles = calloc(sub->n_alloced, sizeof(*sub->particles));
  sub->particles_alt = calloc(sub->n_alloced, sizeof(*sub->particles_alt));
  sub->b_idx = calloc(sub->n_alloced, sizeof(*sub->b_idx));
  sub->b_ids = calloc(sub->n_alloced, sizeof(*sub->b_ids));

  for (int d = 0; d < 3; d++) {
    sub->b_mx[d] = ppsc->patch[prts->p].ldims[d];
    sub->b_dxi[d] = 1.f / ppsc->patch[prts->p].dx[d];
  }
  sub->nr_blocks = sub->b_mx[0] * sub->b_mx[1] * sub->b_mx[2];
  sub->b_cnt = calloc(sub->nr_blocks + 1, sizeof(*sub->b_cnt));
}

static void
psc_particles_single_by_block_destroy(struct psc_particles *prts)
{
  struct psc_particles_single_by_block *sub = psc_particles_single_by_block(prts);

  free(sub->particles);
  free(sub->particles_alt);
  free(sub->b_idx);
  free(sub->b_ids);
  free(sub->b_cnt);
}

static void
psc_particles_single_by_block_reorder(struct psc_particles *prts)
{
  struct psc_particles_single_by_block *sub = psc_particles_single_by_block(prts);

  if (!sub->need_reorder) {
    return;
  }

  for (int n = 0; n < prts->n_part; n++) {
    sub->particles_alt[n] = sub->particles[sub->b_ids[n]];
  }
  
  // swap in alt array
  particle_single_by_block_t *tmp = sub->particles;
  sub->particles = sub->particles_alt;
  sub->particles_alt = tmp;
  sub->need_reorder = false;
}

void
particles_single_by_block_realloc(struct psc_particles *prts, int new_n_part)
{
  struct psc_particles_single_by_block *sub = psc_particles_single_by_block(prts);

  if (new_n_part <= sub->n_alloced)
    return;

  sub->n_alloced = new_n_part * 1.2;
  sub->particles = realloc(sub->particles, sub->n_alloced * sizeof(*sub->particles));
  sub->b_idx = realloc(sub->b_idx, sub->n_alloced * sizeof(*sub->b_idx));
  sub->b_ids = realloc(sub->b_ids, sub->n_alloced * sizeof(*sub->b_ids));
  free(sub->particles_alt);
  sub->particles_alt = malloc(sub->n_alloced * sizeof(*sub->particles_alt));
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
// psc_particles_single_by_block_write

static void
psc_particles_single_by_block_write(struct psc_particles *prts, struct mrc_io *io)
{
  int ierr;
  assert(sizeof(particle_single_by_block_t) / sizeof(particle_single_by_block_real_t) == 8);
  assert(sizeof(particle_single_by_block_real_t) == sizeof(float));

  long h5_file;
  mrc_io_get_h5_file(io, &h5_file);

  hid_t group = H5Gopen(h5_file, mrc_io_obj_path(io, prts), H5P_DEFAULT); H5_CHK(group);
  // save/restore n_alloced, too?
  ierr = H5LTset_attribute_int(group, ".", "p", &prts->p, 1); CE;
  ierr = H5LTset_attribute_int(group, ".", "n_part", &prts->n_part, 1); CE;
  ierr = H5LTset_attribute_uint(group, ".", "flags", &prts->flags, 1); CE;
  if (prts->n_part > 0) {
    // in a rather ugly way, we write the int "kind" member as a float
    hsize_t hdims[2] = { prts->n_part, 8 };
    ierr = H5LTmake_dataset_float(group, "particles_single_by_block", 2, hdims,
				  (float *) particles_single_by_block_get_one(prts, 0)); CE;
  }
  ierr = H5Gclose(group); CE;
}

// ----------------------------------------------------------------------
// psc_particles_single_by_block_by_block_read

static void
psc_particles_single_by_block_read(struct psc_particles *prts, struct mrc_io *io)
{
  int ierr;
  long h5_file;
  mrc_io_get_h5_file(io, &h5_file);

  hid_t group = H5Gopen(h5_file, mrc_io_obj_path(io, prts), H5P_DEFAULT); H5_CHK(group);
  ierr = H5LTget_attribute_int(group, ".", "p", &prts->p); CE;
  ierr = H5LTget_attribute_int(group, ".", "n_part", &prts->n_part); CE;
  ierr = H5LTget_attribute_uint(group, ".", "flags", &prts->flags); CE;
  psc_particles_setup(prts);
  if (prts->n_part > 0) {
    ierr = H5LTread_dataset_float(group, "particles_single_by_block",
				  (float *) particles_single_by_block_get_one(prts, 0)); CE;
  }
  ierr = H5Gclose(group); CE;
}

#endif

// ======================================================================

static void
psc_particles_single_by_block_copy_to_single(struct psc_particles *prts,
					     struct psc_particles *prts_single,
					     unsigned int flags)
{
  prts_single->n_part = prts->n_part;
  assert(prts_single->n_part <= psc_particles_single(prts_single)->n_alloced);
  for (int n = 0; n < prts->n_part; n++) {
    particle_single_by_block_t *prt = particles_single_by_block_get_one(prts, n);
    particle_single_t *prt_single = particles_single_get_one(prts_single, n);
    
    prt_single->xi      = prt->xi;
    prt_single->yi      = prt->yi;
    prt_single->zi      = prt->zi;
    prt_single->pxi     = prt->pxi;
    prt_single->pyi     = prt->pyi;
    prt_single->pzi     = prt->pzi;
    prt_single->qni_wni = prt->qni_wni;
    prt_single->kind    = prt->kind;
  }
}

static void
psc_particles_single_by_block_copy_from_single(struct psc_particles *prts,
					       struct psc_particles *prts_single,
					       unsigned int flags)
{
  prts->n_part = prts_single->n_part;
  assert(prts->n_part <= psc_particles_single_by_block(prts)->n_alloced);
  for (int n = 0; n < prts->n_part; n++) {
    particle_single_by_block_t *prt = particles_single_by_block_get_one(prts, n);
    particle_single_t *prt_single = particles_single_get_one(prts_single, n);
    
    prt->xi      = prt_single->xi;
    prt->yi      = prt_single->yi;
    prt->zi      = prt_single->zi;
    prt->pxi     = prt_single->pxi;
    prt->pyi     = prt_single->pyi;
    prt->pzi     = prt_single->pzi;
    prt->qni_wni = prt_single->qni_wni;
    prt->kind    = prt_single->kind;
  }
}

// ======================================================================
// psc_mparticles: subclass "single_by_block"
  
struct psc_mparticles_ops psc_mparticles_single_by_block_ops = {
  .name                    = "single_by_block",
};

// ======================================================================
// psc_particles: subclass "single"

static struct mrc_obj_method psc_particles_single_by_block_methods[] = {
  MRC_OBJ_METHOD("copy_to_single"  , psc_particles_single_by_block_copy_to_single),
  MRC_OBJ_METHOD("copy_from_single", psc_particles_single_by_block_copy_from_single),
  {}
};

struct psc_particles_ops psc_particles_single_by_block_ops = {
  .name                    = "single_by_block",
  .size                    = sizeof(struct psc_particles_single_by_block),
  .methods                 = psc_particles_single_by_block_methods,
  .setup                   = psc_particles_single_by_block_setup,
  .destroy                 = psc_particles_single_by_block_destroy,
#ifdef HAVE_LIBHDF5_HL
  .read                    = psc_particles_single_by_block_read,
  .write                   = psc_particles_single_by_block_write,
#endif
  .reorder                 = psc_particles_single_by_block_reorder,
};
