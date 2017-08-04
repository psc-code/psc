
#include "psc_particles_single_by_block.h"
#include "psc_particles_as_single_by_block.h"

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

  prts->n_alloced = psc_particles_size(prts) * 1.2;
  sub->particles = calloc(prts->n_alloced, sizeof(*sub->particles));
  sub->particles_alt = calloc(prts->n_alloced, sizeof(*sub->particles_alt));
  sub->b_idx = calloc(prts->n_alloced, sizeof(*sub->b_idx));
  sub->b_ids = calloc(prts->n_alloced, sizeof(*sub->b_ids));

  for (int d = 0; d < 3; d++) {
    sub->b_mx[d] = ppsc->patch[prts->p].ldims[d];
    sub->b_dxi[d] = 1.f / ppsc->patch[prts->p].dx[d];
  }
  sub->nr_blocks = sub->b_mx[0] * sub->b_mx[1] * sub->b_mx[2];
  sub->b_cnt = calloc(sub->nr_blocks + 1, sizeof(*sub->b_cnt));
  sub->b_off = calloc(sub->nr_blocks + 2, sizeof(*sub->b_off));
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
  free(sub->b_off);
}

static void
psc_particles_single_by_block_reorder(struct psc_particles *prts)
{
  struct psc_particles_single_by_block *sub = psc_particles_single_by_block(prts);

  if (!sub->need_reorder) {
    return;
  }

  for (int n = 0; n < psc_particles_size(prts); n++) {
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

  if (new_n_part <= prts->n_alloced)
    return;

  prts->n_alloced = new_n_part * 1.2;
  sub->particles = realloc(sub->particles, prts->n_alloced * sizeof(*sub->particles));
  sub->b_idx = realloc(sub->b_idx, prts->n_alloced * sizeof(*sub->b_idx));
  sub->b_ids = realloc(sub->b_ids, prts->n_alloced * sizeof(*sub->b_ids));
  free(sub->particles_alt);
  sub->particles_alt = malloc(prts->n_alloced * sizeof(*sub->particles_alt));
}

// ----------------------------------------------------------------------
// find_idx_off_1st_rel
//
// FIXME, duplicated and only need fixed shift, no og here

static inline void
find_idx_off_1st_rel(particle_real_t xi[3], int lg[3], particle_real_t og[3], particle_real_t shift,
		     particle_real_t dxi[3])
{
  for (int d = 0; d < 3; d++) {
    particle_real_t pos = xi[d] * dxi[d] + shift;
    lg[d] = particle_real_fint(pos);
    og[d] = pos - lg[d];
  }
}

// ----------------------------------------------------------------------
// psc_particles_single_by_block_get_b_idx

static unsigned int
psc_particles_single_by_block_get_b_idx(struct psc_particles *prts, int n)
{
  struct psc_particles_single_by_block *sub = psc_particles_single_by_block(prts);

  particle_single_by_block_t *prt = particles_single_by_block_get_one(prts, n);
  particle_real_t of[3];
  int b_pos[3], *b_mx = sub->b_mx;
  find_idx_off_1st_rel(&prt->xi, b_pos, of, 0.f, sub->b_dxi);

  // FIXME, only if blocksize == 1, 3D!!!
  if (b_pos[0] >= 0 && b_pos[0] < b_mx[0] &&
      b_pos[1] >= 0 && b_pos[1] < b_mx[1] &&
      b_pos[2] >= 0 && b_pos[2] < b_mx[2]) {
    unsigned int b_idx = (b_pos[2] * b_mx[1] + b_pos[1]) * b_mx[0] + b_pos[0];
    assert(b_idx < sub->nr_blocks);
    return b_idx;
  } else { // out of bounds
    return sub->nr_blocks;
  }
}

// ----------------------------------------------------------------------
// psc_particles_single_by_block_check
//
// make sure that all particles are within the local patch (ie.,
// they have a valid block idx, and that they are sorted by
// block_idx, and that that b_off properly contains the range of
// particles in each block.

static void
psc_particles_single_by_block_check(struct psc_particles *prts)
{
  struct psc_particles_single_by_block *sub = psc_particles_single_by_block(prts);

  assert(psc_particles_size(prts) <= prts->n_alloced);

  int block = 0;
  for (int n = 0; n < psc_particles_size(prts); n++) {
    while (n >= sub->b_off[block + 1]) {
      block++;
      assert(block < sub->nr_blocks);
    }
    assert(n >= sub->b_off[block] && n < sub->b_off[block + 1]);
    assert(block < sub->nr_blocks);
    unsigned int b_idx = psc_particles_single_by_block_get_b_idx(prts, n);
    assert(b_idx < sub->nr_blocks);
    assert(b_idx == block);
  }
}

// ----------------------------------------------------------------------
// psc_particles_single_by_block_sort
//
// This function will take an unsorted list of particles and
// sort it by block index, and set up the b_off array that
// describes the particle index range for each block.
//
// b_idx, b_ids and b_cnt are not valid after returning
// (they shouldn't be needed for anything, either)

static void
psc_particles_single_by_block_sort(struct psc_particles *prts)
{
  struct psc_particles_single_by_block *sub = psc_particles_single_by_block(prts);

  for (int b = 0; b < sub->nr_blocks; b++) {
    sub->b_cnt[b] = 0;
  }

  // calculate block indices for each particle and count
  int n_prts = psc_particles_size(prts);
  for (int n = 0; n < n_prts; n++) {
    unsigned int b_idx = psc_particles_single_by_block_get_b_idx(prts, n);
    assert(b_idx < sub->nr_blocks);
    sub->b_idx[n] = b_idx;
    sub->b_cnt[b_idx]++;
  }
  
  // b_off = prefix sum(b_cnt), zero b_cnt
  int sum = 0;
  for (int b = 0; b <= sub->nr_blocks; b++) {
    sub->b_off[b] = sum;
    sum += sub->b_cnt[b];
    // temporarily abuse b_cnt for next step
    // (the array will be altered, so we use b_cnt rather than
    // b_off, which we'd like to preserve)
    sub->b_cnt[b] = sub->b_off[b];
  }

  // find target position for each particle
  for (int n = 0; n < n_prts; n++) {
    unsigned int b_idx = sub->b_idx[n];
    sub->b_ids[n] = sub->b_cnt[b_idx]++;
  }

  // reorder into alt particle array
  // WARNING: This is reversed to what reorder() does!
  for (int n = 0; n < n_prts; n++) {
    sub->particles_alt[sub->b_ids[n]] = sub->particles[n];
  }
  
  // swap in alt array
  particle_single_by_block_t *tmp = sub->particles;
  sub->particles = sub->particles_alt;
  sub->particles_alt = tmp;
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
  int n_prts = psc_particles_size(prts);
  ierr = H5LTset_attribute_int(group, ".", "n_part", &n_prts, 1); CE;
  ierr = H5LTset_attribute_uint(group, ".", "flags", &prts->flags, 1); CE;
  if (n_prts > 0) {
    // in a rather ugly way, we write the int "kind" member as a float
    hsize_t hdims[2] = { n_prts, 8 };
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
  int n_prts;
  ierr = H5LTget_attribute_int(group, ".", "n_part", &n_prts); CE;
  psc_particles_resize(prts, n_prts);
  ierr = H5LTget_attribute_uint(group, ".", "flags", &prts->flags); CE;
  psc_particles_setup(prts);
  if (n_prts > 0) {
    ierr = H5LTread_dataset_float(group, "particles_single_by_block",
				  (float *) particles_single_by_block_get_one(prts, 0)); CE;
  }
  ierr = H5Gclose(group); CE;
}

#endif

// ======================================================================

static void
psc_mparticles_single_by_block_copy_to_single(int p, struct psc_mparticles *mprts,
					      struct psc_mparticles *mprts_single,
					      unsigned int flags)
{
  struct psc_particles *prts_single = psc_mparticles_get_patch(mprts_single, p);
  struct psc_particles *prts = psc_mparticles_get_patch(mprts, p);
  int n_prts = psc_particles_size(prts);
  psc_particles_resize(prts_single, n_prts);
  assert(n_prts <= prts_single->n_alloced);
  for (int n = 0; n < n_prts; n++) {
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
psc_mparticles_single_by_block_copy_from_single(int p, struct psc_mparticles *mprts,
						struct psc_mparticles *mprts_single,
						unsigned int flags)
{
  struct psc_particles *prts_single = psc_mparticles_get_patch(mprts_single, p);
  struct psc_particles *prts = psc_mparticles_get_patch(mprts, p);
  int n_prts = psc_particles_size(prts_single);
  psc_particles_resize(prts, n_prts);
  assert(n_prts <= prts->n_alloced);
  for (int n = 0; n < n_prts; n++) {
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

  psc_particles_single_by_block_sort(prts);
  psc_particles_single_by_block_check(prts);
}

// ======================================================================
// psc_particles: subclass "single"

struct psc_particles_ops psc_particles_single_by_block_ops = {
  .name                    = "single_by_block",
  .size                    = sizeof(struct psc_particles_single_by_block),
  .setup                   = psc_particles_single_by_block_setup,
  .destroy                 = psc_particles_single_by_block_destroy,
#ifdef HAVE_LIBHDF5_HL
  .read                    = psc_particles_single_by_block_read,
  .write                   = psc_particles_single_by_block_write,
#endif
  .reorder                 = psc_particles_single_by_block_reorder,
};

// ======================================================================
// psc_mparticles: subclass "single_by_block"
  
static struct mrc_obj_method psc_particles_single_by_block_methods[] = {
  MRC_OBJ_METHOD("copy_to_single"  , psc_mparticles_single_by_block_copy_to_single),
  MRC_OBJ_METHOD("copy_from_single", psc_mparticles_single_by_block_copy_from_single),
  {}
};

struct psc_mparticles_ops psc_mparticles_single_by_block_ops = {
  .name                    = "single_by_block",
  .methods                 = psc_particles_single_by_block_methods,
};

