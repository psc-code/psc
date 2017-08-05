
#include "psc.h"
#include "psc_particles_as_single_by_block.h"
#include "psc_particles_inc.h"
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
// conversion to/from "single"

static void
put_particle_single(particle_single_by_block_t *prt, int n, struct psc_particles *prts_dbl)
{
  particle_single_t *prt_dbl = particles_single_get_one(prts_dbl, n);
  
  prt_dbl->xi      = prt->xi;
  prt_dbl->yi      = prt->yi;
  prt_dbl->zi      = prt->zi;
  prt_dbl->pxi     = prt->pxi;
  prt_dbl->pyi     = prt->pyi;
  prt_dbl->pzi     = prt->pzi;
  prt_dbl->qni_wni = prt->qni_wni;
  prt_dbl->kind    = prt->kind;
}

static void
get_particle_single(particle_single_by_block_t *prt, int n, struct psc_particles *prts_dbl)
{
  particle_single_t *prt_dbl = particles_single_get_one(prts_dbl, n);

  prt->xi      = prt_dbl->xi;
  prt->yi      = prt_dbl->yi;
  prt->zi      = prt_dbl->zi;
  prt->pxi     = prt_dbl->pxi;
  prt->pyi     = prt_dbl->pyi;
  prt->pzi     = prt_dbl->pzi;
  prt->qni_wni = prt_dbl->qni_wni;
  prt->kind    = prt_dbl->kind;
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
  psc_particles_single_by_block_sort(psc_mparticles_get_patch(mprts, p));
  psc_particles_single_by_block_check(psc_mparticles_get_patch(mprts, p));
}

// ======================================================================
// psc_particles: subclass "single"

struct psc_particles_ops psc_particles_single_by_block_ops = {
  .name                    = "single_by_block",
  .size                    = sizeof(struct psc_particles_single_by_block),
  .setup                   = psc_particles_single_by_block_setup,
  .destroy                 = psc_particles_single_by_block_destroy,
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

