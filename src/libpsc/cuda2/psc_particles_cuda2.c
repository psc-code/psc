
#include "psc.h"
#include "psc_particles_cuda2.h"
#include "psc_cuda2.h"

// for conversions
#include "psc_particles_single.h"
#include "psc_particles_cuda.h"
#include "../cuda/psc_cuda.h"

#include <stdlib.h>

// ======================================================================
// psc_particles "cuda2"

// ----------------------------------------------------------------------
// psc_particles_cuda2_setup

static void
psc_particles_cuda2_setup(struct psc_particles *prts)
{
  struct psc_particles_cuda2 *sub = psc_particles_cuda2(prts);

  sub->n_alloced = prts->n_part * 1.2;
  sub->particles = calloc(sub->n_alloced, sizeof(*sub->particles));
  sub->particles_alt = calloc(sub->n_alloced, sizeof(*sub->particles_alt));
  sub->b_idx = calloc(sub->n_alloced, sizeof(*sub->b_idx));
  sub->b_ids = calloc(sub->n_alloced, sizeof(*sub->b_ids));

  for (int d = 0; d < 3; d++) {
    sub->dxi[d] = 1.f / ppsc->patch[prts->p].dx[d];
    assert(ppsc->patch[prts->p].ldims[d] % psc_particles_cuda2_bs[d] == 0);
    sub->b_mx[d] = ppsc->patch[prts->p].ldims[d] / psc_particles_cuda2_bs[d];
  }
  sub->nr_blocks = sub->b_mx[0] * sub->b_mx[1] * sub->b_mx[2];
  sub->b_cnt = calloc(sub->nr_blocks + 1, sizeof(*sub->b_cnt));
  sub->b_off = calloc(sub->nr_blocks + 2, sizeof(*sub->b_off));
}

// ----------------------------------------------------------------------
// psc_particles_cuda2_destroy

static void
psc_particles_cuda2_destroy(struct psc_particles *prts)
{
  struct psc_particles_cuda2 *sub = psc_particles_cuda2(prts);

  free(sub->particles);
  free(sub->particles_alt);
  free(sub->b_idx);
  free(sub->b_ids);
  free(sub->b_cnt);
  free(sub->b_off);
}

// ----------------------------------------------------------------------
// find_idx_off_1st_rel
//
// FIXME, duplicated and only need fixed shift, no og here

static inline void
find_idx_off_1st_rel(particle_cuda2_real_t xi[3], int lg[3], particle_cuda2_real_t og[3],
		     particle_cuda2_real_t shift, particle_cuda2_real_t dxi[3])
{
  for (int d = 0; d < 3; d++) {
    particle_cuda2_real_t pos = xi[d] * dxi[d] + shift;
    lg[d] = particle_cuda2_real_fint(pos);
    og[d] = pos - lg[d];
  }
}

// ----------------------------------------------------------------------
// psc_particles_cuda2_get_b_idx

static unsigned int
psc_particles_cuda2_get_b_idx(struct psc_particles *prts, int n)
{
  // FIXME, a particle which ends up exactly on the high boundary
  // should not be considered oob -- well, it's a complicated business,
  // but at least for certain b.c.s it's legal
  struct psc_particles_cuda2 *sub = psc_particles_cuda2(prts);

  particle_cuda2_t *prt = particles_cuda2_get_one(prts, n);
  particle_cuda2_real_t of[3];
  int b_pos[3], *b_mx = sub->b_mx;
  find_idx_off_1st_rel(&prt->xi, b_pos, of, 0.f, sub->dxi);
  for (int d = 0; d < 3; d++) {
    b_pos[d] /= psc_particles_cuda2_bs[d];
  }

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
// psc_particles_cuda2_sort
//
// This function will take an unsorted list of particles and
// sort it by block index, and set up the b_off array that
// describes the particle index range for each block.
//
// b_idx, b_ids and b_cnt are not valid after returning
// (they shouldn't be needed for anything, either)

static void
psc_particles_cuda2_sort(struct psc_particles *prts)
{
  struct psc_particles_cuda2 *sub = psc_particles_cuda2(prts);

  for (int b = 0; b < sub->nr_blocks; b++) {
    sub->b_cnt[b] = 0;
  }

  // calculate block indices for each particle and count
  for (int n = 0; n < prts->n_part; n++) {
    unsigned int b_idx = psc_particles_cuda2_get_b_idx(prts, n);
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
  sub->b_off[sub->nr_blocks + 1] = sum;

  // find target position for each particle
  for (int n = 0; n < prts->n_part; n++) {
    unsigned int b_idx = sub->b_idx[n];
    sub->b_ids[n] = sub->b_cnt[b_idx]++;
  }

  // reorder into alt particle array
  // WARNING: This is reversed to what reorder() does!
  for (int n = 0; n < prts->n_part; n++) {
    sub->particles_alt[sub->b_ids[n]] = sub->particles[n];
  }
  
  // swap in alt array
  particle_cuda2_t *tmp = sub->particles;
  sub->particles = sub->particles_alt;
  sub->particles_alt = tmp;
}

// ----------------------------------------------------------------------
// psc_particles_cuda2_check
//
// make sure that all particles are within the local patch (ie.,
// they have a valid block idx, and that they are sorted by
// block_idx, and that that b_off properly contains the range of
// particles in each block.

static void
psc_particles_cuda2_check(struct psc_particles *prts)
{
  struct psc_particles_cuda2 *sub = psc_particles_cuda2(prts);

  assert(prts->n_part <= sub->n_alloced);

  int block = 0;
  for (int n = 0; n < prts->n_part; n++) {
    while (n >= sub->b_off[block + 1]) {
      block++;
      assert(block < sub->nr_blocks);
    }
    assert(n >= sub->b_off[block] && n < sub->b_off[block + 1]);
    assert(block < sub->nr_blocks);
    unsigned int b_idx = psc_particles_cuda2_get_b_idx(prts, n);
    assert(b_idx < sub->nr_blocks);
    assert(b_idx == block);
  }
}

// ----------------------------------------------------------------------
// psc_particles_cuda2_copy_to_single

static void
psc_particles_cuda2_copy_to_single(struct psc_particles *prts_base,
				    struct psc_particles *prts, unsigned int flags)
{
  prts->n_part = prts_base->n_part;
  assert(prts->n_part <= psc_particles_single(prts)->n_alloced);
  for (int n = 0; n < prts_base->n_part; n++) {
    particle_cuda2_t *part_base = particles_cuda2_get_one(prts_base, n);
    particle_single_t *part = particles_single_get_one(prts, n);
    
    part->xi      = part_base->xi;
    part->yi      = part_base->yi;
    part->zi      = part_base->zi;
    part->pxi     = part_base->pxi;
    part->pyi     = part_base->pyi;
    part->pzi     = part_base->pzi;
    part->qni_wni = part_base->qni_wni;
    part->kind    = part_base->kind;
  }
}

// ----------------------------------------------------------------------
// psc_particles_cuda2_copy_from_single

static void
psc_particles_cuda2_copy_from_single(struct psc_particles *prts_base,
				     struct psc_particles *prts, unsigned int flags)
{
  struct psc_particles_cuda2 *sub = psc_particles_cuda2(prts_base);
  prts_base->n_part = prts->n_part;
  assert(prts_base->n_part <= sub->n_alloced);
  for (int n = 0; n < prts_base->n_part; n++) {
    particle_cuda2_t *part_base = particles_cuda2_get_one(prts_base, n);
    particle_single_t *part = particles_single_get_one(prts, n);

    part_base->xi      = part->xi;
    part_base->yi      = part->yi;
    part_base->zi      = part->zi;
    part_base->pxi     = part->pxi;
    part_base->pyi     = part->pyi;
    part_base->pzi     = part->pzi;
    part_base->qni_wni = part->qni_wni;
    part_base->kind    = part->kind;
  }

  psc_particles_cuda2_sort(prts_base);
  // make sure there are no oob particles
  assert(sub->b_off[sub->nr_blocks] == sub->b_off[sub->nr_blocks+1]);
  psc_particles_cuda2_check(prts_base);
}

// ----------------------------------------------------------------------
// psc_particles_cuda2_copy_to_cuda

static void
psc_particles_cuda2_copy_to_cuda(struct psc_particles *prts,
				 struct psc_particles *prts_cuda, unsigned int flags)
{
  assert(prts_cuda->n_part == prts->n_part);
  
  float4 *xi4  = calloc(prts->n_part, sizeof(float4));
  float4 *pxi4 = calloc(prts->n_part, sizeof(float4));
  
  for (int n = 0; n < prts->n_part; n++) {
    particle_cuda2_t *prt = particles_cuda2_get_one(prts, n);
    
    xi4[n].x  = prt->xi;
    xi4[n].y  = prt->yi;
    xi4[n].z  = prt->zi;
    xi4[n].w  = cuda_int_as_float(prt->kind);
    pxi4[n].x = prt->pxi;
    pxi4[n].y = prt->pyi;
    pxi4[n].z = prt->pzi;
    pxi4[n].w = prt->qni_wni;
  }
  
  __particles_cuda_to_device(prts_cuda, xi4, pxi4);
  
  free(xi4);
  free(pxi4);
}

static void
psc_particles_cuda2_copy_from_cuda(struct psc_particles *prts,
				   struct psc_particles *prts_cuda, unsigned int flags)
{
  struct psc_particles_cuda2 *sub = psc_particles_cuda2(prts);
  prts->n_part = prts_cuda->n_part;
  assert(prts->n_part <= sub->n_alloced);
  
  float4 *xi4  = calloc(prts_cuda->n_part, sizeof(float4));
  float4 *pxi4 = calloc(prts_cuda->n_part, sizeof(float4));
  
  __particles_cuda_from_device(prts_cuda, xi4, pxi4);
  
  for (int n = 0; n < prts->n_part; n++) {
    particle_cuda2_t *prt = particles_cuda2_get_one(prts, n);

    prt->xi      = xi4[n].x;
    prt->yi      = xi4[n].y;
    prt->zi      = xi4[n].z;
    prt->kind    = cuda_float_as_int(xi4[n].w);
    prt->pxi     = pxi4[n].x;
    prt->pyi     = pxi4[n].y;
    prt->pzi     = pxi4[n].z;
    prt->qni_wni = pxi4[n].w;
  }

  free(xi4);
  free(pxi4);
}

// ======================================================================
// psc_particles: subclass "cuda2"

static struct mrc_obj_method psc_particles_cuda2_methods[] = {
  MRC_OBJ_METHOD("copy_to_single"  , psc_particles_cuda2_copy_to_single),
  MRC_OBJ_METHOD("copy_from_single", psc_particles_cuda2_copy_from_single),
  MRC_OBJ_METHOD("copy_to_cuda"    , psc_particles_cuda2_copy_to_cuda),
  MRC_OBJ_METHOD("copy_from_cuda"  , psc_particles_cuda2_copy_from_cuda),
  {}
};

struct psc_particles_ops psc_particles_cuda2_ops = {
  .name                    = "cuda2",
  .size                    = sizeof(struct psc_particles_cuda2),
  .methods                 = psc_particles_cuda2_methods,
  .setup                   = psc_particles_cuda2_setup,
  .destroy                 = psc_particles_cuda2_destroy,
#if 0
#ifdef HAVE_LIBHDF5_HL
  .read                    = psc_particles_cuda2_read,
  .write                   = psc_particles_cuda2_write,
#endif
  .reorder                 = psc_particles_cuda2_reorder,
#endif
};

// ======================================================================
// psc_mparticles: subclass "cuda2"
  
struct psc_mparticles_ops psc_mparticles_cuda2_ops = {
  .name                    = "cuda2",
};

