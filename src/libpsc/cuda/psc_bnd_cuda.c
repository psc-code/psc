
#include "psc_bnd_cuda.h"
#include "psc_cuda.h"
#include "psc_bnd_cuda_fields.h"
#include "psc_particles_as_single.h"
#include "../psc_bnd/ddc_particles.h"

#include <mrc_profile.h>

struct psc_bnd_sub {
};

#define to_psc_bnd_sub(bnd) ((struct psc_bnd_sub *)((bnd)->obj.subctx))

struct cuda_ctx_patch {
  particle_t *prts;
  float4 *xi4;
  float4 *pxi4;
  int n;
};

struct cuda_ctx {
  struct cuda_ctx_patch *cuda_patch;
};

// ----------------------------------------------------------------------
// ddcp_particles helpers

static void
ddcp_particles_realloc(void *_ctx, int p, int new_n_particles)
{
  struct cuda_ctx *ctx = _ctx;

  assert(!ctx->cuda_patch[p].prts);
  ctx->cuda_patch[p].prts = malloc(new_n_particles * sizeof(*ctx->cuda_patch[p].prts));
}

static void *
ddcp_particles_get_addr(void *_ctx, int p, int n)
{
  struct cuda_ctx *ctx = _ctx;
  return ctx->cuda_patch[p].prts + n;
}

// ----------------------------------------------------------------------
// psc_bnd_sub_setup

static void
psc_bnd_sub_setup(struct psc_bnd *bnd)
{
  psc_bnd_setup_super(bnd);
  bnd->ddcp = ddc_particles_create(bnd->ddc, sizeof(particle_t),
				   sizeof(particle_real_t),
				   MPI_PARTICLES_REAL,
				   ddcp_particles_realloc,
				   ddcp_particles_get_addr);
}

// ----------------------------------------------------------------------
// psc_bnd_sub_unsetup

static void
psc_bnd_sub_unsetup(struct psc_bnd *bnd)
{
  ddc_particles_destroy(bnd->ddcp);
}

// ----------------------------------------------------------------------
// cpatch_append helper

static void
cpatch_append(struct psc_particles *prts, struct cuda_ctx_patch *cpatch,
	      particle_t *prt,
	      unsigned int *bn_idx, unsigned int *bn_cnts, unsigned int *bn_off)
{
  struct psc_particles_cuda *cuda = psc_particles_cuda(prts);
  int nn = cpatch->n++;
  cpatch->xi4[nn].x  = prt->xi;
  cpatch->xi4[nn].y  = prt->yi;
  cpatch->xi4[nn].z  = prt->zi;
  cpatch->xi4[nn].w  = cuda_int_as_float(prt->kind);
  cpatch->pxi4[nn].x = prt->pxi;
  cpatch->pxi4[nn].y = prt->pyi;
  cpatch->pxi4[nn].z = prt->pzi;
  cpatch->pxi4[nn].w = prt->qni_wni;

  int b_pos[3];
  for (int d = 0; d < 3; d++) {
    float *xi = &cpatch->xi4[nn].x;
    b_pos[d] = particle_real_fint(xi[d] * cuda->b_dxi[d]);
    if (b_pos[d] < 0 || b_pos[d] >= cuda->b_mx[d]) {
      printf("!!! xi %g %g %g\n", xi[0], xi[1], xi[2]);
      printf("!!! d %d xi4[n] %g biy %d // %d\n",
	     d, xi[d], b_pos[d], cuda->b_mx[d]);
      if (b_pos[d] < 0) {
	xi[d] = 0.f;
      } else {
	xi[d] *= (1. - 1e-6);
      }
    }
    b_pos[d] = particle_real_fint(xi[d] * cuda->b_dxi[d]);
    assert(b_pos[d] >= 0 && b_pos[d] < cuda->b_mx[d]);
  }
  if (bn_idx) {
    unsigned int b =
      (b_pos[2] * cuda->b_mx[1] + b_pos[1]) * cuda->b_mx[0] + b_pos[0];
    assert(b < cuda->nr_blocks);
    bn_idx[nn] = b;
    bn_off[nn] = bn_cnts[b]++;
  }
}

// ----------------------------------------------------------------------
// exchange_particles_host
//
// Does the CPU part of the particle communication,
// starting at the copy from the GPU, exchange,
// copy back to GPU.
// Calculates counts and offsets for the newly received particles
// to help the GPU put them into the right final place.

static void
exchange_particles_host(struct psc_bnd *bnd, struct psc_mparticles *particles,
			int *offsets, int *n_sends,
			unsigned int **bn_cnts, unsigned int **bn_idx,
			unsigned int **bn_off)
{
  struct psc *psc = bnd->psc;
  struct ddc_particles *ddcp = bnd->ddcp;

  static int pr_A, pr_B, pr_D;
  if (!pr_A) {
    pr_A = prof_register("xchg_prep_" PARTICLE_TYPE, 1., 0, 0);
    pr_B = prof_register("xchg_comm_" PARTICLE_TYPE, 1., 0, 0);
    pr_D = prof_register("xchg_to_dev", 1., 0, 0);
  }

  prof_start(pr_A);
  struct cuda_ctx ctx;
  ctx.cuda_patch = calloc(particles->nr_patches, sizeof(*ctx.cuda_patch));
  for (int p = 0; p < particles->nr_patches; p++) {
    struct psc_particles *prts_cuda = psc_mparticles_get_patch(particles, p);

    struct cuda_ctx_patch *cpatch = &ctx.cuda_patch[p];
    cpatch->n = 0;
    cpatch->prts = NULL;

    int n_send = n_sends[p];
    cpatch->xi4  = malloc(n_send * sizeof(*cpatch->xi4));
    cpatch->pxi4 = malloc(n_send * sizeof(*cpatch->pxi4));
    if (bn_idx) {
      bn_idx[p] = malloc(n_send * sizeof(*bn_idx[p]));
      bn_off[p] = malloc(n_send * sizeof(*bn_off[p]));
    }

    __particles_cuda_from_device_range(prts_cuda, cpatch->xi4, cpatch->pxi4,
				       offsets[p], offsets[p] + n_send);
  }

  for (int p = 0; p < particles->nr_patches; p++) {
    struct psc_particles *prts_cuda = psc_mparticles_get_patch(particles, p);
    struct psc_particles_cuda *cuda = psc_particles_cuda(prts_cuda);
    struct cuda_ctx_patch *cpatch = &ctx.cuda_patch[p];
    struct ddcp_patch *ddcp_patch = &ddcp->patches[p];
    ddcp_patch->head = 0;
    for (int dir1 = 0; dir1 < N_DIR; dir1++) {
      ddcp_patch->nei[dir1].n_send = 0;
    }
    struct psc_patch *patch = psc->patch + p;
    particle_real_t xm[3];
    for (int d = 0; d < 3; d++) {
      xm[d] = patch->ldims[d] * psc->dx[d];
    }

    particle_real_t *b_dxi = cuda->b_dxi;
    int *b_mx = cuda->b_mx;

    for (int n = 0; n < n_sends[p]; n++) {
      particle_t prt;
      prt.xi      = cpatch->xi4[n].x;
      prt.yi      = cpatch->xi4[n].y;
      prt.zi      = cpatch->xi4[n].z;
      prt.kind    = cuda_float_as_int(cpatch->xi4[n].w);
      prt.pxi     = cpatch->pxi4[n].x;
      prt.pyi     = cpatch->pxi4[n].y;
      prt.pzi     = cpatch->pxi4[n].z;
      prt.qni_wni = cpatch->pxi4[n].w;

      particle_real_t *xi = &prt.xi;
      particle_real_t *pxi = &prt.pxi;

      bool drop = false;
      int dir[3];
      for (int d = 0; d < 3; d++) {
	int bi = particle_real_fint(xi[d] * b_dxi[d]);
	if (bi < 0) {
	  // FIXME, assumes every patch has same dimensions
	  if (patch->off[d] != 0 || psc->domain.bnd_part_lo[d] == BND_PART_PERIODIC) {
	    xi[d] += xm[d];
	    dir[d] = -1;
	    bi = particle_real_fint(xi[d] * b_dxi[d]);
	    if (bi >= b_mx[d]) {
	      xi[d] = 0.;
	      dir[d] = 0;
	    }
	  } else {
	    switch (psc->domain.bnd_part_lo[d]) {
	    case BND_PART_REFLECTING:
	      xi[d]  = -xi[d];
	      pxi[d] = -pxi[d];
	      dir[d] = 0;
	      break;
	    case BND_PART_ABSORBING:
	      drop = true;
	      break;
	    default:
	      assert(0);
	    }
	  }
	} else if (bi >= b_mx[d]) {
	  if (patch->off[d] + patch->ldims[d] != psc->domain.gdims[d] ||
	      psc->domain.bnd_part_hi[d] == BND_PART_PERIODIC) {
	    xi[d] -= xm[d];
	    dir[d] = +1;
	    bi = particle_real_fint(xi[d] * b_dxi[d]);
	    if (bi < 0) {
	      xi[d] = 0.;
	    }
	  } else {
	    switch (psc->domain.bnd_part_hi[d]) {
	    case BND_PART_REFLECTING:
	      xi[d] = 2.f * xm[d] - xi[d];
	      pxi[d] = -pxi[d];
	      dir[d] = 0;
	      bi = particle_real_fint(xi[d] * b_dxi[d]);
	      if (bi >= b_mx[d]) {
		xi[d] *= (1. - 1e-6);
	      }
	      break;
	    case BND_PART_ABSORBING:
	      drop = true;
	      break;
	    default:
	      assert(0);
	    }
	  }
	} else {
	  dir[d] = 0;
	}
	assert(xi[d] >= 0.f);
	assert(xi[d] <= xm[d]);
      }
      if (!drop) {
	if (dir[0] == 0 && dir[1] == 0 && dir[2] == 0) {
	  cpatch_append(prts_cuda, cpatch, &prt, bn_idx[p], bn_cnts[p], bn_off[p]);
	} else {
	  ddc_particles_queue(ddcp, ddcp_patch, dir, &prt);
	}
      }
    }
  }
  prof_stop(pr_A);

  prof_start(pr_B);
  ddc_particles_comm(ddcp, &ctx);
  prof_stop(pr_B);

  prof_start(pr_D);
  for (int p = 0; p < particles->nr_patches; p++) {
    struct ddcp_patch *patch = &ddcp->patches[p];
    struct cuda_ctx_patch *cpatch = ctx.cuda_patch + p;
    struct psc_particles *prts_cuda = psc_mparticles_get_patch(particles, p);
    struct psc_particles_cuda *cuda = psc_particles_cuda(prts_cuda);
    int n_recv = cpatch->n + patch->head;
    prts_cuda->n_part = offsets[p] + n_recv;
    assert(prts_cuda->n_part <= cuda->n_alloced);

    cpatch->xi4  = realloc(cpatch->xi4, n_recv * sizeof(float4));
    cpatch->pxi4 = realloc(cpatch->pxi4, n_recv * sizeof(float4));
    bn_idx[p] = realloc(bn_idx[p], n_recv * sizeof(*bn_idx[p]));
    bn_off[p] = realloc(bn_off[p], n_recv * sizeof(*bn_off[p]));
    for (int n = 0; n < patch->head; n++) {
      cpatch_append(prts_cuda, cpatch, &cpatch->prts[n], bn_idx[p], bn_cnts[p], bn_off[p]);
    }

    __particles_cuda_to_device_range(prts_cuda, cpatch->xi4, cpatch->pxi4,
				     offsets[p], offsets[p] + n_recv);

    free(cpatch->prts);
    free(cpatch->xi4);
    free(cpatch->pxi4);
  }
  prof_stop(pr_D);

  free(ctx.cuda_patch);
}

// ----------------------------------------------------------------------
// psc_bnd_sub_exchange_particles_serial_periodic
//
// specialized version if there's only one single patch,
// with all periodic b.c.
//
// TODO: a lot of consolidation with the generic case should
// be possible. In particular, optimizations carry over, like
// calculating new offsets as part of the spine, not calculating
// new block indices, not calculates old ids.
// The most significant stumbling block to make them pretty much the same
// is that this case here needs to handle the odd periodic spine.

static void
psc_bnd_sub_exchange_particles_serial_periodic(struct psc_bnd *psc_bnd,
						mparticles_cuda_t *particles)
{
  static int pr, pr_F, pr_G, pr_H;
  if (!pr) {
    pr   = prof_register("xchg_parts", 1., 0, 0);
    pr_F = prof_register("xchg_bidx_ids", 1., 0, 0);
    pr_G = prof_register("xchg_sort_pairs", 1., 0, 0);
    pr_H = prof_register("xchg_reorder_off", 1., 0, 0);
  }

  prof_start(pr);
  cuda_exchange_particles(0, psc_mparticles_get_patch(particles, 0));

  // sort
  for (int p = 0; p < particles->nr_patches; p++) {
    struct psc_particles *prts = psc_mparticles_get_patch(particles, p);
    struct psc_particles_cuda *cuda = psc_particles_cuda(prts);

    prof_start(pr_F);
    cuda_find_block_indices(prts, cuda->d_part.bidx);
    prof_stop(pr_F);

    prof_start(pr_G);
    sort_pairs_device_2(cuda->d_part.sort_ctx, cuda->d_part.bidx,
			cuda->d_part.alt_ids,
			prts->n_part,
			cuda->d_part.offsets);
    prof_stop(pr_G);

    prof_start(pr_H);
    cuda_reorder(prts, cuda->d_part.alt_ids);
    prof_stop(pr_H);
  }
  prof_stop(pr);
}

// ----------------------------------------------------------------------
// psc_bnd_sub_exchange_particles_general

static void
psc_bnd_sub_exchange_particles_general(struct psc_bnd *psc_bnd,
					mparticles_cuda_t *particles)
{
  static int pr, pr_B, pr_C, pr_D, pr_F, pr_G, pr_H;
  if (!pr) {
    pr   = prof_register("xchg_parts", 1., 0, 0);
    pr_B = prof_register("xchg_bidx", 1., 0, 0);
    pr_C = prof_register("xchg_pfxsum", 1., 0, 0);
    pr_D = prof_register("xchg_reorder", 1., 0, 0);
    pr_F = prof_register("xchg_bidx_ids", 1., 0, 0);
    pr_G = prof_register("xchg_sort_pairs", 1., 0, 0);
    pr_H = prof_register("xchg_reorder_off", 1., 0, 0);
  }
  
  prof_start(pr);
  int *n_sends = malloc(ppsc->nr_patches * sizeof(*n_sends));
  int *n_parts = malloc(ppsc->nr_patches * sizeof(*n_parts));
  unsigned int **bn_cnts = malloc(ppsc->nr_patches * sizeof(*bn_cnts));
  unsigned int **bn_idx = malloc(ppsc->nr_patches * sizeof(*bn_idx));
  unsigned int **bn_off = malloc(ppsc->nr_patches * sizeof(*bn_off));
  for (int p = 0; p < particles->nr_patches; p++) {
    struct psc_particles *prts = psc_mparticles_get_patch(particles, p);
    struct psc_particles_cuda *cuda = psc_particles_cuda(prts);
    bn_cnts[p] = calloc(cuda->nr_blocks, sizeof(*bn_cnts[p]));

    prof_start(pr_B);
    cuda_find_block_indices_2(prts, cuda->d_part.bidx, 0);
    prof_stop(pr_B);

    prof_start(pr_C);
    n_sends[p] = cuda_exclusive_scan_2(prts, cuda->d_part.bidx, cuda->d_part.sums);
    n_parts[p] = prts->n_part;
    prof_stop(pr_C);

    prof_start(pr_D);
    cuda_reorder_send_buf(p, prts, cuda->d_part.bidx, cuda->d_part.sums, n_sends[p]);
    prof_stop(pr_D);
  }

  exchange_particles_host(psc_bnd, particles, n_parts, n_sends,
			  bn_cnts, bn_idx, bn_off);

  for (int p = 0; p < particles->nr_patches; p++) {
    struct psc_particles *prts = psc_mparticles_get_patch(particles, p);
    struct psc_particles_cuda *cuda = psc_particles_cuda(prts);

    prof_start(pr_F);
    cuda_find_block_indices_3(prts, cuda->d_part.bidx, cuda->d_part.alt_bidx,
			      n_parts[p], bn_idx[p], bn_off[p]);
    free(bn_idx[p]);
    free(bn_off[p]);
    prof_stop(pr_F);

    // OPT: when calculating bidx, do preprocess then
    prof_start(pr_G);
    void *sp = sort_pairs_3_create(cuda->b_mx);
    sort_pairs_3_device(sp, cuda->d_part.bidx, cuda->d_part.alt_bidx, cuda->d_part.alt_ids,
			prts->n_part, cuda->d_part.offsets,
			n_parts[p], bn_cnts[p]);
    sort_pairs_3_destroy(sp);
    prof_stop(pr_G);

    prof_start(pr_H);
    cuda_reorder(prts, cuda->d_part.alt_ids);
    prof_stop(pr_H);

    prts->n_part = n_parts[p] - n_sends[p] + prts->n_part - n_parts[p];

    free(bn_cnts[p]);
  }

  free(n_sends);
  free(n_parts);
  free(bn_cnts);
  free(bn_idx);
  free(bn_off);
  prof_stop(pr);
}

// ----------------------------------------------------------------------
// psc_bnd_sub_exchange_particles

void
psc_bnd_sub_exchange_particles(struct psc_bnd *bnd,
				mparticles_base_t *particles_base)
{
  int size;
  MPI_Comm_size(psc_bnd_comm(bnd), &size);

  // This function only makes sense if it's called for particles already being of cuda
  // type. We could call _get_cuda(), but that wouldn't be happy if some particles were
  // not in the right patch in the first place.

  assert(strcmp(psc_mparticles_type(particles_base), "cuda") == 0);
  mparticles_cuda_t *particles = particles_base;

  if (size == 1 && ppsc->nr_patches == 1 && // FIXME !!!
      ppsc->domain.bnd_fld_lo[0] == BND_FLD_PERIODIC &&
      ppsc->domain.bnd_fld_lo[1] == BND_FLD_PERIODIC &&
      ppsc->domain.bnd_fld_lo[2] == BND_FLD_PERIODIC) {
    psc_bnd_sub_exchange_particles_serial_periodic(bnd, particles);
  } else {
    psc_bnd_sub_exchange_particles_general(bnd, particles);
  }
}

// ======================================================================
// psc_bnd: subclass "cuda"

struct psc_bnd_ops psc_bnd_cuda_ops = {
  .name                  = "cuda",
  .size                  = sizeof(struct psc_bnd_sub),
  .setup                 = psc_bnd_sub_setup,
  .unsetup               = psc_bnd_sub_unsetup,
  .exchange_particles    = psc_bnd_sub_exchange_particles,

  .create_ddc            = psc_bnd_fields_cuda_create,
  .add_ghosts            = psc_bnd_fields_cuda_add_ghosts,
  .fill_ghosts           = psc_bnd_fields_cuda_fill_ghosts,
};

