
#include "psc.h"
#include "psc_particles_acc.h"
#include "psc_acc.h"

// for conversions
#include "psc_particles_single.h"

#include <stdlib.h>

// ======================================================================
// psc_particles "acc"

// ----------------------------------------------------------------------
// find_idx_off_1st_rel
//
// FIXME, duplicated and only need fixed shift, no og here

static inline void
find_idx_off_1st_rel(particle_acc_real_t xi[3], int lg[3], particle_acc_real_t og[3],
		     particle_acc_real_t shift, particle_acc_real_t dxi[3])
{
  for (int d = 0; d < 3; d++) {
    particle_acc_real_t pos = xi[d] * dxi[d] + shift;
    lg[d] = fint(pos);
    og[d] = pos - lg[d];
  }
}

// ----------------------------------------------------------------------
// psc_particles_acc_get_b_idx

static unsigned int
psc_particles_acc_get_b_idx(struct psc_particles *prts, int n)
{
  // FIXME, a particle which ends up exactly on the high boundary
  // should not be considered oob -- well, it's a complicated business,
  // but at least for certain b.c.s it's legal
  struct psc_particles_acc *sub = psc_particles_acc(prts);

  particle_acc_t prt;
  PARTICLE_ACC_LOAD_POS(prt, sub->xi4, n);
  PARTICLE_ACC_LOAD_MOM(prt, sub->pxi4, n);
  particle_acc_real_t of[3];
  int b_pos[3], *b_mx = sub->b_mx;
  find_idx_off_1st_rel(&particle_acc_x(&prt), b_pos, of, 0.f, sub->dxi);
  for (int d = 0; d < 3; d++) {
    b_pos[d] /= psc_particles_acc_bs[d];
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
// psc_particles_acc_sort
//
// This function will take an unsorted list of particles and
// sort it by block index, and set up the b_off array that
// describes the particle index range for each block.
//
// b_idx, b_ids and b_cnt are not valid after returning
// (they shouldn't be needed for anything, either)

static void
psc_particles_acc_sort(struct psc_particles *prts)
{
  struct psc_particles_acc *sub = psc_particles_acc(prts);

  for (int b = 0; b < sub->nr_blocks; b++) {
    sub->b_cnt[b] = 0;
  }

  // calculate block indices for each particle and count
  for (int n = 0; n < prts->n_part; n++) {
    unsigned int b_idx = psc_particles_acc_get_b_idx(prts, n);
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
    sub->xi4_alt [sub->b_ids[n]] = sub->xi4 [n];
    sub->pxi4_alt[sub->b_ids[n]] = sub->pxi4[n];
  }
  
  // swap in alt array
  float4 *tmp_xi4 = sub->xi4;
  float4 *tmp_pxi4 = sub->pxi4;
  sub->xi4 = sub->xi4_alt;
  sub->pxi4 = sub->pxi4_alt;
  sub->xi4_alt = tmp_xi4;
  sub->pxi4_alt = tmp_pxi4;
}

// ----------------------------------------------------------------------
// psc_particles_acc_check
//
// make sure that all particles are within the local patch (ie.,
// they have a valid block idx, and that they are sorted by
// block_idx, and that that b_off properly contains the range of
// particles in each block.

static void
psc_particles_acc_check(struct psc_particles *prts)
{
  struct psc_particles_acc *sub = psc_particles_acc(prts);

  assert(prts->n_part <= sub->n_alloced);

  int block = 0;
  for (int n = 0; n < prts->n_part; n++) {
    while (n >= sub->b_off[block + 1]) {
      block++;
      assert(block < sub->nr_blocks);
    }
    assert(n >= sub->b_off[block] && n < sub->b_off[block + 1]);
    assert(block < sub->nr_blocks);
    unsigned int b_idx = psc_particles_acc_get_b_idx(prts, n);
    assert(b_idx < sub->nr_blocks);
    assert(b_idx == block);
  }
}

static void
copy_from(int p, struct psc_mparticles *mprts,
	  struct psc_mparticles *mprts_acc, unsigned int flags,
	  void (*get_particle)(particle_double_t *prt, int n, struct psc_particles *prts))
{
  struct psc_particles *prts = psc_mparticles_get_patch(mprts, p);
  struct psc_particles *prts_acc = psc_mparticles_get_patch(mprts_acc, p);
  int n_prts = psc_particles_size(prts_acc);
  mprts[p].resize(n_prts);
  for (int n = 0; n < n_prts; n++) {
    particle_acc_t prt;
    get_particle(&prt, n, prts_acc);
    PARTICLE_ACC_STORE_POS(prt, sub->xi4, n);
    PARTICLE_ACC_STORE_MOM(prt, sub->pxi4, n);
  }
}

static void
copy_to(int p, struct psc_mparticles *mprts,
	struct psc_mparticles *mprts_acc, unsigned int flags,
	void (*put_particle)(particle_double_t *prt, int n, struct psc_particles *prts))
{
  struct psc_particles *prts = psc_mparticles_get_patch(mprts, p);
  struct psc_particles *prts_acc = psc_mparticles_get_patch(mprts_acc, p);
  int n_prts = psc_particles_size(prts);
  mprts_acc[p].resize(n_prts);
  for (int n = 0; n < n_prts; n++) {
    particle_acc_t prt;
    PARTICLE_ACC_LOAD_POS(prt, sub->xi4, n);
    PARTICLE_ACC_LOAD_MOM(prt, sub->pxi4, n);
    put_particle(&prt, n, prts_acc);
  }
}

// ======================================================================
// conversion to/from "single"

static void
get_particle_single(particle_acc_t *prt, int n, struct psc_particles *prts_single)
{
  particle_single_t *prt_single = particles_single_get_one(prts_single, n);

  prt->xi      = prt_single->xi;
  prt->yi      = prt_single->yi;
  prt->zi      = prt_single->zi;
  prt->kind    = prt_single->kind;
  prt->pxi     = prt_single->pxi;
  prt->pyi     = prt_single->pyi;
  prt->pzi     = prt_single->pzi;
  prt->qni_wni = prt_single->qni_wni;
}

static void
put_particle_single(particle_acc_t *prt, int n, struct psc_particles *prts_single)
{
  particle_single_t *prt_single = particles_single_get_one(prts_single, n);

  prt_single->xi      = prt->xi;
  prt_single->yi      = prt->yi;
  prt_single->zi      = prt->zi;
  prt_single->kind    = prt->kind;
  prt_single->pxi     = prt->pxi;
  prt_single->pyi     = prt->pyi;
  prt_single->pzi     = prt->pzi;
  prt_single->qni_wni = prt->qni_wni;
}

static void
psc_particles_acc_copy_to_single(int p, struct psc_mparticles *mprts_base,
				 struct psc_mparticles *mprts, unsigned int flags)
{
  copy_to(p, mprts_base, mprts, flags, put_particle_single);
}

static void
psc_particles_acc_copy_from_single(int p, struct psc_mparticles *mprts_base,
				   struct psc_mparticles *mprts, unsigned int flags)
{
  copy_from(p, mprts_base, mprts, flags, get_particle_single);
  }
}

// ----------------------------------------------------------------------
// psc_particles: subclass "acc"

static struct mrc_obj_method psc_particles_acc_methods[] = {
  MRC_OBJ_METHOD("copy_to_single"  , psc_particles_acc_copy_to_single),
  MRC_OBJ_METHOD("copy_from_single", psc_particles_acc_copy_from_single),
  {}
};

// ======================================================================
// psc_mparticles: subclass "acc"
  
// ----------------------------------------------------------------------
// psc_mparticles_acc_setup

static void
psc_mparticles_acc_setup(struct psc_mparticles *mprts)
{
  struct psc_mparticles_acc *sub = psc_mparticles_acc(mprts);

  psc_mparticles_setup_super(mprts);

  if (mprts->nr_patches == 0) {
    return;
  }
  
  int *gdims = ppsc->domain.gdims;
  for (int d = 0; d < 3; d++) {
    sub->bs[d] = psc_particles_acc_bs[d];
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
    struct psc_particles *prts = psc_mparticles_get_patch(mprts, p);

    sub->n_part_total += prts->n_part;
  }
  sub->n_alloced_total = sub->n_part_total * 1.2;

  sub->xi4 = calloc(sub->n_alloced_total, sizeof(*sub->xi4));
  sub->pxi4 = calloc(sub->n_alloced_total, sizeof(*sub->pxi4));
  sub->xi4_alt = calloc(sub->n_alloced_total, sizeof(*sub->xi4_alt));
  sub->pxi4_alt = calloc(sub->n_alloced_total, sizeof(*sub->pxi4_alt));
  sub->b_off = calloc(sub->nr_blocks_total + 2, sizeof(*sub->b_off));
  
  float4 *xi4 = sub->xi4;
  float4 *pxi4 = sub->pxi4;
  float4 *xi4_alt = sub->xi4_alt;
  float4 *pxi4_alt = sub->pxi4_alt;
  for (int p = 0; p < mprts->nr_patches; p++) {
    struct psc_particles *prts = psc_mparticles_get_patch(mprts, p);
    struct psc_particles_acc *prts_sub = psc_particles_acc(prts);

    prts_sub->n_alloced = prts->n_part;

    for (int d = 0; d < 3; d++) {
      prts_sub->dxi[d] = sub->dxi[d];
      prts_sub->b_mx[d] = sub->b_mx[d];
    }
    prts_sub->nr_blocks = sub->nr_blocks;

    // on host
    prts_sub->xi4 = xi4;
    prts_sub->pxi4 = pxi4;
    prts_sub->xi4_alt = xi4_alt;
    prts_sub->pxi4_alt = pxi4_alt;
    xi4 += prts->n_part;
    pxi4 += prts->n_part;
    xi4_alt += prts->n_part;
    pxi4_alt += prts->n_part;

    prts_sub->b_idx = calloc(prts_sub->n_alloced, sizeof(*prts_sub->b_idx));
    prts_sub->b_ids = calloc(prts_sub->n_alloced, sizeof(*prts_sub->b_ids));
    prts_sub->b_cnt = calloc(prts_sub->nr_blocks + 1, sizeof(*prts_sub->b_cnt));
    prts_sub->b_off = calloc(prts_sub->nr_blocks + 2, sizeof(*prts_sub->b_off));
  }
}

// ----------------------------------------------------------------------
// psc_mparticles_acc_destroy

static void
psc_mparticles_acc_destroy(struct psc_mparticles *mprts)
{
  struct psc_mparticles_acc *sub = psc_mparticles_acc(mprts);

  for (int p = 0; p < mprts->nr_patches; p++) {
    struct psc_particles *prts = psc_mparticles_get_patch(mprts, p);
    struct psc_particles_acc *prts_sub = psc_particles_acc(prts);

    free(prts_sub->b_idx);
    free(prts_sub->b_ids);
    free(prts_sub->b_cnt);
    free(prts_sub->b_off);
  }

  free(sub->xi4);
  free(sub->pxi4);
  free(sub->xi4_alt);
  free(sub->pxi4_alt);
  free(sub->b_off);
}

// ----------------------------------------------------------------------
// psc_mparticles_acc_check
//
// make sure that h_b_off[] is set up properly

static void
psc_mparticles_acc_check(struct psc_mparticles *mprts)
{
  struct psc_mparticles_acc *sub = psc_mparticles_acc(mprts);

  assert(sub->n_part_total <= sub->n_alloced_total);

  int block = 0;
  for (int n = 0; n < sub->n_part_total; n++) {
    while (n >= sub->b_off[block + 1]) {
      block++;
      assert(block < sub->nr_blocks_total);
    }
    assert(n >= sub->b_off[block] && n < sub->b_off[block + 1]);
    assert(block < sub->nr_blocks_total);
#if 0
    unsigned int b_idx = psc_particles_acc_get_b_idx(prts, n);
    assert(b_idx < sub->nr_blocks);
    assert(b_idx == block);
#endif
  }
}

// ----------------------------------------------------------------------
// psc_mparticles_acc_setup_internals

static void
psc_mparticles_acc_setup_internals(struct psc_mparticles *mprts)
{
  struct psc_mparticles_acc *sub = psc_mparticles_acc(mprts);

  int nr_blocks = sub->nr_blocks;
  // FIXME, should just do the sorting all-in-one
  for (int p = 0; p < mprts->nr_patches; p++) {
    struct psc_particles *prts = psc_mparticles_get_patch(mprts, p);
    struct psc_particles_acc *prts_sub = psc_particles_acc(prts);

    psc_particles_acc_sort(prts);
    // make sure there are no oob particles
    assert(prts_sub->b_off[nr_blocks] == prts_sub->b_off[nr_blocks+1]);
    psc_particles_acc_check(prts);
  }
  // FIXME, to keep consistency with per-patch
  // swap in alt array
  float4 *tmp_xi4 = sub->xi4;
  float4 *tmp_pxi4 = sub->pxi4;
  sub->xi4 = sub->xi4_alt;
  sub->pxi4 = sub->pxi4_alt;
  sub->xi4_alt = tmp_xi4;
  sub->pxi4_alt = tmp_pxi4;

  unsigned int n_part = 0;
  for (int p = 0; p < mprts->nr_patches; p++) {
    struct psc_particles *prts = psc_mparticles_get_patch(mprts, p);
    struct psc_particles_acc *prts_sub = psc_particles_acc(prts);

    for (int b = 0; b < nr_blocks; b++) {
      sub->b_off[p * nr_blocks + b] = prts_sub->b_off[b] + n_part;
    }
    n_part += prts->n_part;
  }
  sub->b_off[sub->nr_blocks_total] = sub->n_part_total;
  sub->b_off[sub->nr_blocks_total + 1] = sub->n_part_total;

  psc_mparticles_acc_check(mprts);
}

// ----------------------------------------------------------------------
// psc_mparticles: subclass "acc"

struct psc_mparticles_ops psc_mparticles_acc_ops = {
  .name                    = "acc",
  .size                    = sizeof(struct psc_mparticles_acc),
  .setup                   = psc_mparticles_acc_setup,
  .destroy                 = psc_mparticles_acc_destroy,
  .setup_internals         = psc_mparticles_acc_setup_internals,
};

