
#include "psc.h"
#include "psc_particles_cuda2.h"
#include "psc_cuda2.h"

// for conversions
#include "psc_particles_single.h"
#include "psc_particles_cuda.h"

#include "../cuda/cuda_mparticles.h"

#include <stdlib.h>

// ======================================================================
// psc_particles "cuda2"

// ----------------------------------------------------------------------
// psc_mparticles_cuda2_copy_to_device

void
psc_mparticles_cuda2_copy_to_device(struct psc_mparticles *mprts)
{
  struct psc_mparticles_cuda2 *sub = psc_mparticles_cuda2(mprts);

  cuda_memcpy_device_from_host(sub->d_xi4, sub->h_xi4, sub->n_part_total * sizeof(*sub->d_xi4));
  cuda_memcpy_device_from_host(sub->d_pxi4, sub->h_pxi4, sub->n_part_total * sizeof(*sub->d_pxi4));
  cuda_memcpy_device_from_host(sub->d_b_off, sub->h_b_off, (sub->nr_blocks_total + 2) * sizeof(*sub->d_b_off));
}

// ----------------------------------------------------------------------
// psc_mparticles_cuda2_copy_to_host

void
psc_mparticles_cuda2_copy_to_host(struct psc_mparticles *mprts)
{
  struct psc_mparticles_cuda2 *sub = psc_mparticles_cuda2(mprts);

  cuda_memcpy_host_from_device(sub->h_xi4, sub->d_xi4, sub->n_part_total * sizeof(*sub->h_xi4));
  cuda_memcpy_host_from_device(sub->h_pxi4, sub->d_pxi4, sub->n_part_total * sizeof(*sub->h_pxi4));
  cuda_memcpy_host_from_device(sub->h_b_off, sub->d_b_off, (sub->nr_blocks_total + 2) * sizeof(*sub->h_b_off));
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
    lg[d] = fint(pos);
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

  particle_cuda2_t prt;
  PARTICLE_CUDA2_LOAD_POS(prt, sub->h_xi4, n);
  PARTICLE_CUDA2_LOAD_MOM(prt, sub->h_pxi4, n);
  particle_cuda2_real_t of[3];
  int b_pos[3], *b_mx = sub->b_mx;
  find_idx_off_1st_rel(&particle_cuda2_x(&prt), b_pos, of, 0.f, sub->dxi);
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
  struct psc_mparticles *mprts = prts->mprts;
  int p = prts->p;
  struct psc_particles_cuda2 *sub = psc_particles_cuda2(prts);

  for (int b = 0; b < sub->nr_blocks; b++) {
    sub->b_cnt[b] = 0;
  }

  // calculate block indices for each particle and count
  int n_prts = psc_mparticles_n_prts_by_patch(mprts, p);
  for (int n = 0; n < n_prts; n++) {
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
  for (int n = 0; n < n_prts; n++) {
    unsigned int b_idx = sub->b_idx[n];
    sub->b_ids[n] = sub->b_cnt[b_idx]++;
  }

  // reorder into alt particle array
  // WARNING: This is reversed to what reorder() does!
  for (int n = 0; n < n_prts; n++) {
    sub->h_xi4_alt [sub->b_ids[n]] = sub->h_xi4 [n];
    sub->h_pxi4_alt[sub->b_ids[n]] = sub->h_pxi4[n];
  }
  
  // swap in alt array
  float4 *tmp_xi4 = sub->h_xi4;
  float4 *tmp_pxi4 = sub->h_pxi4;
  sub->h_xi4 = sub->h_xi4_alt;
  sub->h_pxi4 = sub->h_pxi4_alt;
  sub->h_xi4_alt = tmp_xi4;
  sub->h_pxi4_alt = tmp_pxi4;
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
  struct psc_mparticles *mprts = prts->mprts;
  int p = prts->p;
  struct psc_particles_cuda2 *sub = psc_particles_cuda2(prts);

  int n_prts = psc_mparticles_n_prts_by_patch(mprts, p);
  assert(n_prts <= psc_mparticles_n_alloced(mprts, p));

  int block = 0;
  for (int n = 0; n < n_prts; n++) {
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
  struct psc_mparticles *mprts = prts->mprts;
  int p = prts->p;
  struct psc_particles_cuda2 *sub = psc_particles_cuda2(prts_base);

  int n_prts = psc_mparticles_n_prts_by_patch(prts_base->mprts, p);
  mprts[p].resize(n_prts);
  for (int n = 0; n < n_prts; n++) {
    particle_cuda2_t prt_base;
    PARTICLE_CUDA2_LOAD_POS(prt_base, sub->h_xi4, n);
    PARTICLE_CUDA2_LOAD_MOM(prt_base, sub->h_pxi4, n);
    particle_single_t *part = psc_mparticles_single_get_one(prts->mprts, prts->p, n);
    
    part->xi      = prt_base.xi[0];
    part->yi      = prt_base.xi[1];
    part->zi      = prt_base.xi[2];
    part->kind    = particle_cuda2_kind(&prt_base);
    part->pxi     = prt_base.pxi[0];
    part->pyi     = prt_base.pxi[1];
    part->pzi     = prt_base.pxi[2];
    part->qni_wni = prt_base.qni_wni;
  }
}

// ----------------------------------------------------------------------
// psc_particles_cuda2_copy_from_single

static void
psc_particles_cuda2_copy_from_single(struct psc_particles *prts_base,
				     struct psc_particles *prts, unsigned int flags)
{
  struct psc_particles_cuda2 *sub = psc_particles_cuda2(prts_base);

  int n_prts = psc_mparticles_n_prts_by_patch(prts->mprts, prts->p);
  prts_base->mprts[prts_base->p].resize(n_prts);
  for (int n = 0; n < n_prts; n++) {
    particle_cuda2_t prt_base;
    particle_single_t *part = psc_mparticles_single_get_one(prts->mprts, prts->p, n);

    prt_base.xi[0]         = part->xi;
    prt_base.xi[1]         = part->yi;
    prt_base.xi[2]         = part->zi;
    prt_base.kind_as_float = cuda_int_as_float(part->kind);
    prt_base.pxi[0]        = part->pxi;
    prt_base.pxi[1]        = part->pyi;
    prt_base.pxi[2]        = part->pzi;
    prt_base.qni_wni       = part->qni_wni;

    PARTICLE_CUDA2_STORE_POS(prt_base, sub->h_xi4, n);
    PARTICLE_CUDA2_STORE_MOM(prt_base, sub->h_pxi4, n);
  }
}

#ifdef USE_CUDA

#include "../cuda/psc_cuda.h"

static void
particles_cuda_to_device(struct psc_particles *prts, float4 *xi4, float4 *pxi4)
{
  struct psc_mparticles *mprts = prts->mprts;
  struct psc_mparticles_cuda *mprts_cuda = psc_mparticles_cuda(mprts);
  struct cuda_mparticles *cmprts = mprts_cuda->cmprts;

  unsigned int off = 0;
  for (int p = 0; p < prts->p; p++) {
    off += psc_mparticles_n_prts_by_patch(mprts, p);
  }

  cuda_mparticles_to_device(cmprts, xi4, pxi4, psc_mparticles_n_prts_by_patch(mprts, prts->p), off);
}

// ----------------------------------------------------------------------
// psc_particles_cuda2_copy_to_cuda

static void
psc_particles_cuda2_copy_to_cuda(struct psc_particles *prts,
				 struct psc_particles *prts_cuda, unsigned int flags)
{
  struct psc_particles_cuda2 *sub = psc_particles_cuda2(prts);

  int n_prts = psc_mparticles_n_prts_by_patch(prts->mprts, prts->p);
  assert(psc_mparticles_n_prts_by_patch(prts_cuda->mprts, prts->p) == n_prts);
  
  float4 *xi4  = calloc(n_prts, sizeof(float4));
  float4 *pxi4 = calloc(n_prts, sizeof(float4));
  
  for (int n = 0; n < n_prts; n++) {
    particle_cuda2_t prt;
    PARTICLE_CUDA2_LOAD_POS(prt, sub->h_xi4, n);
    PARTICLE_CUDA2_LOAD_MOM(prt, sub->h_pxi4, n);
    
    xi4[n].x  = prt.xi[0];
    xi4[n].y  = prt.xi[1];
    xi4[n].z  = prt.xi[2];
    xi4[n].w  = prt.kind_as_float;
    pxi4[n].x = prt.pxi[0];
    pxi4[n].y = prt.pxi[1];
    pxi4[n].z = prt.pxi[2];
    pxi4[n].w = prt.qni_wni;
  }
  
  particles_cuda_to_device(prts_cuda, xi4, pxi4);
  
  free(xi4);
  free(pxi4);
}

static void
psc_particles_cuda2_copy_from_cuda(struct psc_particles *prts,
				   struct psc_particles *prts_cuda, unsigned int flags)
{
  struct psc_particles_cuda2 *sub = psc_particles_cuda2(prts);

  int n_prts = psc_mparticles_n_prts_by_patch(prts_cuda->mprts, prts_cuda->p);
  prts->mprts[prts->p].resize(n_prts);
  unsigned int off = 0;
  for (int p = 0; p < prts->p; p++) {
    off += psc_mparticles_n_prts_by_patch(prts_cuda->mprts, p);
  }
  
  float4 *xi4  = calloc(n_prts, sizeof(float4));
  float4 *pxi4 = calloc(n_prts, sizeof(float4));
  
  __particles_cuda_from_device(prts_cuda->mprts, xi4, pxi4, off, n_prts);
  
  for (int n = 0; n < n_prts; n++) {
    particle_cuda2_t prt;

    prt.xi[0]         = xi4[n].x;
    prt.xi[1]         = xi4[n].y;
    prt.xi[2]         = xi4[n].z;
    prt.kind_as_float = xi4[n].w;
    prt.pxi[0]        = pxi4[n].x;
    prt.pxi[1]        = pxi4[n].y;
    prt.pxi[2]        = pxi4[n].z;
    prt.qni_wni       = pxi4[n].w;

    PARTICLE_CUDA2_STORE_POS(prt, sub->h_xi4, n);
    PARTICLE_CUDA2_STORE_MOM(prt, sub->h_pxi4, n);
  }

  free(xi4);
  free(pxi4);
}

#endif

// ----------------------------------------------------------------------
// psc_particles_cuda2_methods

static struct mrc_obj_method psc_particles_cuda2_methods[] = {
  MRC_OBJ_METHOD("copy_to_single"  , psc_particles_cuda2_copy_to_single),
  MRC_OBJ_METHOD("copy_from_single", psc_particles_cuda2_copy_from_single),
#ifdef USE_CUDA
  MRC_OBJ_METHOD("copy_to_cuda"    , psc_particles_cuda2_copy_to_cuda),
  MRC_OBJ_METHOD("copy_from_cuda"  , psc_particles_cuda2_copy_from_cuda),
#endif
  {}
};

// ----------------------------------------------------------------------
// psc_mparticles_cuda2_setup

static void
psc_mparticles_cuda2_setup(struct psc_mparticles *mprts)
{
  struct psc_mparticles_cuda2 *sub = psc_mparticles_cuda2(mprts);

  psc_mparticles_setup_super(mprts);

  if (mprts->nr_patches == 0) {
    return;
  }
  
  int *gdims = ppsc->domain.gdims;
  for (int d = 0; d < 3; d++) {
    sub->bs[d] = psc_particles_cuda2_bs[d];
    if (gdims[d] == 1) {
      sub->bs[d] = 1;
    }
    sub->dxi[d] = 1.f / ppsc->patch[0].dx[d];
    assert(ppsc->patch[0].ldims[d] % sub->bs[d] == 0);
    sub->b_mx[d] = ppsc->patch[0].ldims[d] / sub->bs[d];
  }
  sub->nr_blocks = sub->b_mx[0] * sub->b_mx[1] * sub->b_mx[2];
  sub->nr_blocks_total = sub->nr_blocks * mprts->nr_patches;

  sub->n_part_total = 0;
  for (int p = 0; p < mprts->nr_patches; p++) {
    sub->n_part_total += psc_mparticles_n_prts_by_patch(mprts, p);
  }
  sub->n_alloced_total = sub->n_part_total * 1.2;

  sub->h_xi4 = calloc(sub->n_alloced_total, sizeof(*sub->h_xi4));
  sub->h_pxi4 = calloc(sub->n_alloced_total, sizeof(*sub->h_pxi4));
  sub->h_xi4_alt = calloc(sub->n_alloced_total, sizeof(*sub->h_xi4_alt));
  sub->h_pxi4_alt = calloc(sub->n_alloced_total, sizeof(*sub->h_pxi4_alt));
  sub->h_b_off = calloc(sub->nr_blocks_total + 2, sizeof(*sub->h_b_off));

  sub->d_xi4 = cuda_calloc(sub->n_alloced_total, sizeof(*sub->d_xi4));
  sub->d_pxi4 = cuda_calloc(sub->n_alloced_total, sizeof(*sub->d_pxi4));
  sub->d_b_off = cuda_calloc(sub->nr_blocks_total + 2, sizeof(*sub->d_b_off));
  
  float4 *h_xi4 = sub->h_xi4;
  float4 *h_pxi4 = sub->h_pxi4;
  float4 *h_xi4_alt = sub->h_xi4_alt;
  float4 *h_pxi4_alt = sub->h_pxi4_alt;
  for (int p = 0; p < mprts->nr_patches; p++) {
    struct psc_particles *prts = psc_mparticles_get_patch(mprts, p);
    struct psc_particles_cuda2 *prts_sub = psc_particles_cuda2(prts);

    int n_alloced = psc_mparticles_n_prts_by_patch(mprts, p);
    psc_mparticles_set_n_alloced(mprts, p, n_alloced);

    for (int d = 0; d < 3; d++) {
      prts_sub->dxi[d] = sub->dxi[d];
      prts_sub->b_mx[d] = sub->b_mx[d];
    }
    prts_sub->nr_blocks = sub->nr_blocks;

    // on host
    prts_sub->h_xi4 = h_xi4;
    prts_sub->h_pxi4 = h_pxi4;
    prts_sub->h_xi4_alt = h_xi4_alt;
    prts_sub->h_pxi4_alt = h_pxi4_alt;
    int n_prts = psc_mparticles_n_prts_by_patch(mprts, p);
    h_xi4 += n_prts;
    h_pxi4 += n_prts;
    h_xi4_alt += n_prts;
    h_pxi4_alt += n_prts;

    prts_sub->b_idx = calloc(n_alloced, sizeof(*prts_sub->b_idx));
    prts_sub->b_ids = calloc(n_alloced, sizeof(*prts_sub->b_ids));
    prts_sub->b_cnt = calloc(prts_sub->nr_blocks + 1, sizeof(*prts_sub->b_cnt));
    prts_sub->b_off = calloc(prts_sub->nr_blocks + 2, sizeof(*prts_sub->b_off));
  }
}

// ----------------------------------------------------------------------
// psc_mparticles_cuda2_destroy

static void
psc_mparticles_cuda2_destroy(struct psc_mparticles *mprts)
{
  struct psc_mparticles_cuda2 *sub = psc_mparticles_cuda2(mprts);

  for (int p = 0; p < mprts->nr_patches; p++) {
    struct psc_particles *prts = psc_mparticles_get_patch(mprts, p);
    struct psc_particles_cuda2 *prts_sub = psc_particles_cuda2(prts);

    free(prts_sub->b_idx);
    free(prts_sub->b_ids);
    free(prts_sub->b_cnt);
    free(prts_sub->b_off);
  }

  free(sub->h_xi4);
  free(sub->h_pxi4);
  free(sub->h_xi4_alt);
  free(sub->h_pxi4_alt);
  free(sub->h_b_off);

  cuda_free(sub->d_xi4);
  cuda_free(sub->d_pxi4);
  cuda_free(sub->d_b_off);
}

// ----------------------------------------------------------------------
// psc_mparticles_cuda2_check
//
// make sure that h_b_off[] is set up properly

static void
psc_mparticles_cuda2_check(struct psc_mparticles *mprts)
{
  struct psc_mparticles_cuda2 *sub = psc_mparticles_cuda2(mprts);

  assert(sub->n_part_total <= sub->n_alloced_total);

  int block = 0;
  for (int n = 0; n < sub->n_part_total; n++) {
    while (n >= sub->h_b_off[block + 1]) {
      block++;
      assert(block < sub->nr_blocks_total);
    }
    assert(n >= sub->h_b_off[block] && n < sub->h_b_off[block + 1]);
    assert(block < sub->nr_blocks_total);
#if 0
    unsigned int b_idx = psc_particles_cuda2_get_b_idx(prts, n);
    assert(b_idx < sub->nr_blocks);
    assert(b_idx == block);
#endif
  }
}

// ----------------------------------------------------------------------
// psc_mparticles_cuda2_setup_internals

static void
psc_mparticles_cuda2_setup_internals(struct psc_mparticles *mprts)
{
  struct psc_mparticles_cuda2 *sub = psc_mparticles_cuda2(mprts);

  int nr_blocks = sub->nr_blocks;
  // FIXME, should just do the sorting all-in-one
  for (int p = 0; p < mprts->nr_patches; p++) {
    struct psc_particles *prts = psc_mparticles_get_patch(mprts, p);
    struct psc_particles_cuda2 *prts_sub = psc_particles_cuda2(prts);

    psc_particles_cuda2_sort(prts);
    // make sure there are no oob particles
    assert(prts_sub->b_off[nr_blocks] == prts_sub->b_off[nr_blocks+1]);
    psc_particles_cuda2_check(prts);
  }
  // FIXME, to keep consistency with per-patch
  // swap in alt array
  float4 *tmp_xi4 = sub->h_xi4;
  float4 *tmp_pxi4 = sub->h_pxi4;
  sub->h_xi4 = sub->h_xi4_alt;
  sub->h_pxi4 = sub->h_pxi4_alt;
  sub->h_xi4_alt = tmp_xi4;
  sub->h_pxi4_alt = tmp_pxi4;

  unsigned int n_part = 0;
  for (int p = 0; p < mprts->nr_patches; p++) {
    struct psc_particles *prts = psc_mparticles_get_patch(mprts, p);
    struct psc_particles_cuda2 *prts_sub = psc_particles_cuda2(prts);

    for (int b = 0; b < nr_blocks; b++) {
      sub->h_b_off[p * nr_blocks + b] = prts_sub->b_off[b] + n_part;
    }
    n_part += psc_mparticles_n_prts_by_patch(mprts, p);
  }
  sub->h_b_off[sub->nr_blocks_total] = sub->n_part_total;
  sub->h_b_off[sub->nr_blocks_total + 1] = sub->n_part_total;

  psc_mparticles_cuda2_check(mprts);
}

// ----------------------------------------------------------------------
// psc_mparticles: subclass "cuda2"

struct psc_mparticles_ops psc_mparticles_cuda2_ops = {
  .name                    = "cuda2",
  .size                    = sizeof(struct psc_mparticles_cuda2),
  .methods                 = psc_particles_cuda2_methods,
  .setup                   = psc_mparticles_cuda2_setup,
  .destroy                 = psc_mparticles_cuda2_destroy,
  .setup_internals         = psc_mparticles_cuda2_setup_internals,
};

