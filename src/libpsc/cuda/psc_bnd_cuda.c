
#include "psc_bnd_cuda.h"
#include "psc_cuda.h"
#include "psc_bnd_cuda_fields.h"
#include "psc_particles_as_single.h"
#include "../psc_bnd/ddc_particles.h"

#include <mrc_profile.h>

struct psc_bnd_sub {
};

#define to_psc_bnd_sub(bnd) ((struct psc_bnd_sub *)((bnd)->obj.subctx))

// ----------------------------------------------------------------------
// ddcp_particles helpers

static void
ddcp_particles_realloc(void *_ctx, int p, int new_n_particles)
{
  struct psc_mparticles *mprts = _ctx;
  struct psc_particles *prts = psc_mparticles_get_patch(mprts, p);
  struct psc_particles_cuda *cuda = psc_particles_cuda(prts);

  assert(!cuda->bnd_prts);
  cuda->bnd_prts = malloc(new_n_particles * sizeof(*cuda->bnd_prts));
}

static void *
ddcp_particles_get_addr(void *_ctx, int p, int n)
{
  struct psc_mparticles *mprts = _ctx;
  struct psc_particles *prts = psc_mparticles_get_patch(mprts, p);
  struct psc_particles_cuda *cuda = psc_particles_cuda(prts);

  return &cuda->bnd_prts[n];
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
// xchg_append helper

static void
xchg_append(struct psc_particles *prts, particle_t *prt)
{
  struct psc_particles_cuda *cuda = psc_particles_cuda(prts);
  int nn = cuda->bnd_n_part++;
  cuda->bnd_xi4[nn].x  = prt->xi;
  cuda->bnd_xi4[nn].y  = prt->yi;
  cuda->bnd_xi4[nn].z  = prt->zi;
  cuda->bnd_xi4[nn].w  = cuda_int_as_float(prt->kind);
  cuda->bnd_pxi4[nn].x = prt->pxi;
  cuda->bnd_pxi4[nn].y = prt->pyi;
  cuda->bnd_pxi4[nn].z = prt->pzi;
  cuda->bnd_pxi4[nn].w = prt->qni_wni;

  int b_pos[3];
  for (int d = 0; d < 3; d++) {
    float *xi = &cuda->bnd_xi4[nn].x;
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
  unsigned int b =
    (b_pos[2] * cuda->b_mx[1] + b_pos[1]) * cuda->b_mx[0] + b_pos[0];
  assert(b < cuda->nr_blocks);
  cuda->bnd_idx[nn] = b;
  cuda->bnd_off[nn] = cuda->bnd_cnt[b]++;
}

static inline particle_t *
xchg_get_one(struct psc_particles *prts, int n)
{
  struct psc_particles_cuda *cuda = psc_particles_cuda(prts);
  static particle_t prt_static;

  particle_t *prt = &prt_static;
  prt->xi      = cuda->bnd_xi4[n].x;
  prt->yi      = cuda->bnd_xi4[n].y;
  prt->zi      = cuda->bnd_xi4[n].z;
  prt->kind    = cuda_float_as_int(cuda->bnd_xi4[n].w);
  prt->pxi     = cuda->bnd_pxi4[n].x;
  prt->pyi     = cuda->bnd_pxi4[n].y;
  prt->pzi     = cuda->bnd_pxi4[n].z;
  prt->qni_wni = cuda->bnd_pxi4[n].w;

  return prt;
}

static inline int *
get_b_mx(struct psc_particles *prts)
{
  struct psc_particles_cuda *cuda = psc_particles_cuda(prts);
  return cuda->b_mx;
}

static inline particle_real_t *
get_b_dxi(struct psc_particles *prts)
{
  struct psc_particles_cuda *cuda = psc_particles_cuda(prts);
  return cuda->b_dxi;
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
exchange_particles_host(struct psc_bnd *bnd, struct psc_mparticles *particles)
{
  struct psc *psc = bnd->psc;
  struct ddc_particles *ddcp = bnd->ddcp;

  static int pr_A, pr_B, pr_D;
  if (!pr_A) {
    pr_A = prof_register("xchg_prep_cuda", 1., 0, 0);
    pr_B = prof_register("xchg_comm_cuda", 1., 0, 0);
    pr_D = prof_register("xchg_post", 1., 0, 0);
  }

  prof_start(pr_A);

  // FIXME we should make sure (assert) we don't quietly drop particle which left
  // in the invariant direction

  for (int p = 0; p < particles->nr_patches; p++) {
    struct psc_particles *prts = psc_mparticles_get_patch(particles, p);
    struct psc_particles_cuda *cuda = psc_particles_cuda(prts);

    struct psc_patch *patch = &psc->patch[p];
    particle_real_t xm[3];
    for (int d = 0; d < 3; d++) {
      xm[d] = patch->ldims[d] * psc->dx[d];
    }
    particle_real_t *b_dxi = get_b_dxi(prts);
    int *b_mx = get_b_mx(prts);

    struct ddcp_patch *ddcp_patch = &ddcp->patches[p];
    ddcp_patch->head = 0;
    for (int dir1 = 0; dir1 < N_DIR; dir1++) {
      ddcp_patch->nei[dir1].n_send = 0;
    }
    for (int n = ddcp_patch->head; n < cuda->bnd_n_send; n++) {
      particle_t *prt = xchg_get_one(prts, n);
      particle_real_t *xi = &prt->xi;
      particle_real_t *pxi = &prt->pxi;
      
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
	  xchg_append(prts, prt);
	} else {
	  ddc_particles_queue(ddcp, ddcp_patch, dir, prt);
	}
      }
    }
  }
  prof_stop(pr_A);

  prof_start(pr_B);
  ddc_particles_comm(ddcp, particles);
  prof_stop(pr_B);

  prof_start(pr_D);
  for (int p = 0; p < particles->nr_patches; p++) {
    struct ddcp_patch *patch = &ddcp->patches[p];
    struct psc_particles *prts_cuda = psc_mparticles_get_patch(particles, p);
    struct psc_particles_cuda *cuda = psc_particles_cuda(prts_cuda);
    int n_recv = cuda->bnd_n_part + patch->head;
    prts_cuda->n_part = cuda->bnd_n_part_save + n_recv;
    assert(prts_cuda->n_part <= cuda->n_alloced);

    cuda->bnd_xi4  = realloc(cuda->bnd_xi4, n_recv * sizeof(float4));
    cuda->bnd_pxi4 = realloc(cuda->bnd_pxi4, n_recv * sizeof(float4));
    cuda->bnd_idx  = realloc(cuda->bnd_idx, n_recv * sizeof(*cuda->bnd_idx));
    cuda->bnd_off  = realloc(cuda->bnd_off, n_recv * sizeof(*cuda->bnd_off));
    for (int n = 0; n < patch->head; n++) {
      xchg_append(prts_cuda, &cuda->bnd_prts[n]);
    }

    __particles_cuda_to_device_range(prts_cuda, cuda->bnd_xi4, cuda->bnd_pxi4,
				     cuda->bnd_n_part_save, cuda->bnd_n_part_save + n_recv);

    free(cuda->bnd_prts);
    free(cuda->bnd_xi4);
    free(cuda->bnd_pxi4);
  }
  prof_stop(pr_D);
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
// xchg_copy_from_dev

static void
xchg_copy_from_dev(struct psc_bnd *bnd, struct psc_particles *prts)
{
  struct psc_particles_cuda *cuda = psc_particles_cuda(prts);
  
  cuda->bnd_n_part = 0;
  cuda->bnd_prts = NULL;
  
  int n_send = cuda->bnd_n_send;
  cuda->bnd_xi4  = malloc(n_send * sizeof(*cuda->bnd_xi4));
  cuda->bnd_pxi4 = malloc(n_send * sizeof(*cuda->bnd_pxi4));
  cuda->bnd_idx  = malloc(n_send * sizeof(*cuda->bnd_idx));
  cuda->bnd_off  = malloc(n_send * sizeof(*cuda->bnd_off));
  
  // OPT, could use streaming
  __particles_cuda_from_device_range(prts, cuda->bnd_xi4, cuda->bnd_pxi4,
				     cuda->bnd_n_part_save, cuda->bnd_n_part_save + n_send);
}

// ----------------------------------------------------------------------
// psc_bnd_sub_exchange_particles_general

static void
psc_bnd_sub_exchange_particles_general(struct psc_bnd *psc_bnd,
					mparticles_cuda_t *particles)
{
  static int pr, pr_B, pr_C, pr_D, pr_E, pr_F, pr_G, pr_H;
  if (!pr) {
    pr   = prof_register("xchg_parts", 1., 0, 0);
    pr_B = prof_register("xchg_bidx", 1., 0, 0);
    pr_C = prof_register("xchg_pfxsum", 1., 0, 0);
    pr_D = prof_register("xchg_reorder", 1., 0, 0);
    pr_E = prof_register("xchg_copy_from_dev", 1., 0, 0);
    pr_F = prof_register("xchg_bidx_ids", 1., 0, 0);
    pr_G = prof_register("xchg_sort_pairs", 1., 0, 0);
    pr_H = prof_register("xchg_reorder_off", 1., 0, 0);
  }
  
  prof_start(pr);
  for (int p = 0; p < particles->nr_patches; p++) {
    struct psc_particles *prts = psc_mparticles_get_patch(particles, p);
    struct psc_particles_cuda *cuda = psc_particles_cuda(prts);
    cuda->bnd_cnt = calloc(cuda->nr_blocks, sizeof(*cuda->bnd_cnt));
    
    prof_start(pr_B);
    cuda_find_block_indices_2(prts, cuda->d_part.bidx, 0);
    prof_stop(pr_B);

    prof_start(pr_C);
    cuda->bnd_n_send = cuda_exclusive_scan_2(prts, cuda->d_part.bidx, cuda->d_part.sums);
    cuda->bnd_n_part_save = prts->n_part;
    prof_stop(pr_C);

    prof_start(pr_D);
    cuda_reorder_send_buf(p, prts, cuda->d_part.bidx, cuda->d_part.sums, cuda->bnd_n_send);
    prof_stop(pr_D);

    prof_start(pr_E);
    xchg_copy_from_dev(psc_bnd, prts);
    prof_stop(pr_E);
  }

  exchange_particles_host(psc_bnd, particles);

  for (int p = 0; p < particles->nr_patches; p++) {
    struct psc_particles *prts = psc_mparticles_get_patch(particles, p);
    struct psc_particles_cuda *cuda = psc_particles_cuda(prts);

    prof_start(pr_F);
    cuda_find_block_indices_3(prts, cuda->d_part.bidx, cuda->d_part.alt_bidx,
			      cuda->bnd_n_part_save, cuda->bnd_idx, cuda->bnd_off);
    free(cuda->bnd_idx);
    free(cuda->bnd_off);
    prof_stop(pr_F);

    // OPT: when calculating bidx, do preprocess then
    prof_start(pr_G);
    void *sp = sort_pairs_3_create(cuda->b_mx);
    sort_pairs_3_device(sp, cuda->d_part.bidx, cuda->d_part.alt_bidx, cuda->d_part.alt_ids,
			prts->n_part, cuda->d_part.offsets,
			cuda->bnd_n_part_save, cuda->bnd_cnt);
    sort_pairs_3_destroy(sp);
    prof_stop(pr_G);

    prof_start(pr_H);
    cuda_reorder(prts, cuda->d_part.alt_ids);
    prof_stop(pr_H);

    prts->n_part -= cuda->bnd_n_send;

    free(cuda->bnd_cnt);
  }

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

