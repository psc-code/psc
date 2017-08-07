
#include "psc.h"
#include "psc_particles_as_single_by_block.h"
#include "psc_particles_inc.h"
#include "psc_particles_single.h"

#if 0

static void _mrc_unused // FIXME
psc_particles_single_by_block_reorder(struct psc_mparticles *mprts, int p)
{
  struct psc_mparticles_single_by_block *sub = psc_mparticles_single_by_block(mprts);

  if (!patch->need_reorder) {
    return;
  }

  for (int n = 0; n < psc_mparticles_n_prts_by_patch(mprts, p); n++) {
    patch->prt_array_alt[n] = patch->prt_array[patch->b_ids[n]];
  }
  
  // swap in alt array
  particle_single_by_block_t *tmp = patch->particles;
  patch->particles = patch->particles_alt;
  patch->particles_alt = tmp;
  patch->need_reorder = false;
}

#endif

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
psc_particles_single_by_block_get_b_idx(struct psc_mparticles *mprts, int p, int n)
{
  struct psc_mparticles_single_by_block *msub = psc_mparticles_single_by_block(mprts);
  struct psc_mparticles_single_by_block_patch *patch = &msub->patch[p];
  
  particle_single_by_block_t *prt = psc_mparticles_single_by_block_get_one(mprts, p, n);
  particle_real_t of[3];
  int b_pos[3], *b_mx = patch->b_mx;
  find_idx_off_1st_rel(&prt->xi, b_pos, of, 0.f, patch->b_dxi);

  // FIXME, only if blocksize == 1, 3D!!!
  if (b_pos[0] >= 0 && b_pos[0] < b_mx[0] &&
      b_pos[1] >= 0 && b_pos[1] < b_mx[1] &&
      b_pos[2] >= 0 && b_pos[2] < b_mx[2]) {
    unsigned int b_idx = (b_pos[2] * b_mx[1] + b_pos[1]) * b_mx[0] + b_pos[0];
    assert(b_idx < patch->nr_blocks);
    return b_idx;
  } else { // out of bounds
    return patch->nr_blocks;
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
psc_particles_single_by_block_check(struct psc_mparticles *mprts, int p)
{
  struct psc_mparticles_single_by_block *msub = psc_mparticles_single_by_block(mprts);
  struct psc_mparticles_single_by_block_patch *patch = &msub->patch[p];

  int n_prts = psc_mparticles_n_prts_by_patch(mprts, p);
  assert(n_prts <= patch->n_alloced);

  int block = 0;
  for (int n = 0; n < n_prts; n++) {
    while (n >= patch->b_off[block + 1]) {
      block++;
      assert(block < patch->nr_blocks);
    }
    assert(n >= patch->b_off[block] && n < patch->b_off[block + 1]);
    assert(block < patch->nr_blocks);
    unsigned int b_idx = psc_particles_single_by_block_get_b_idx(mprts, p, n);
    assert(b_idx < patch->nr_blocks);
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
psc_particles_single_by_block_sort(struct psc_mparticles *mprts, int p)
{
  assert(0);
#if 0
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
#endif
}

// ======================================================================
// psc_mparticles: subclass "single_by_block"
  
// ----------------------------------------------------------------------
// conversion to/from "single"

static void
put_particle_single(particle_single_by_block_t *prt, int n, struct psc_mparticles *mprts_sngl, int p)
{
  particle_single_t *prt_sngl = psc_mparticles_single_get_one(mprts_sngl, p, n);
  
  prt_sngl->xi      = prt->xi;
  prt_sngl->yi      = prt->yi;
  prt_sngl->zi      = prt->zi;
  prt_sngl->pxi     = prt->pxi;
  prt_sngl->pyi     = prt->pyi;
  prt_sngl->pzi     = prt->pzi;
  prt_sngl->qni_wni = prt->qni_wni;
  prt_sngl->kind    = prt->kind;
}

static void
get_particle_single(particle_single_by_block_t *prt, int n, struct psc_mparticles *mprts_sngl, int p)
{
  particle_single_t *prt_sngl = psc_mparticles_single_get_one(mprts_sngl, p, n);

  prt->xi      = prt_sngl->xi;
  prt->yi      = prt_sngl->yi;
  prt->zi      = prt_sngl->zi;
  prt->pxi     = prt_sngl->pxi;
  prt->pyi     = prt_sngl->pyi;
  prt->pzi     = prt_sngl->pzi;
  prt->qni_wni = prt_sngl->qni_wni;
  prt->kind    = prt_sngl->kind;
}

static void
psc_mparticles_single_by_block_copy_to_single(int p, struct psc_mparticles *mprts,
				    struct psc_mparticles *mprts_dbl, unsigned int flags)
{
  psc_mparticles_copy_to(p, mprts, mprts_dbl, flags, put_particle_single);
}

static void
psc_mparticles_single_by_block_copy_from_single(int p, struct psc_mparticles *mprts,
				       struct psc_mparticles *mprts_dbl, unsigned int flags)
{
  psc_mparticles_copy_from(p, mprts, mprts_dbl, flags, get_particle_single);
  psc_particles_single_by_block_sort(mprts, p);
  psc_particles_single_by_block_check(mprts, p);
}

static struct mrc_obj_method psc_mparticles_single_by_block_methods[] = {
  MRC_OBJ_METHOD("copy_to_single"  , psc_mparticles_single_by_block_copy_to_single),
  MRC_OBJ_METHOD("copy_from_single", psc_mparticles_single_by_block_copy_from_single),
  {}
};

#include "psc_particles_common.c"

