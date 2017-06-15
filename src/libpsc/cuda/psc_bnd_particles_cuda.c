
#include "psc_bnd_particles_private.h"

#include "psc_bnd_cuda.h"
#include "psc_cuda.h"
#include "particles_cuda.h"
#include "psc_particles_as_single.h"
#include "../psc_bnd/ddc_particles.h"
#include "cuda_mparticles.h"

#include <mrc_profile.h>

// ----------------------------------------------------------------------
// ddcp_particles helpers

static void
ddcp_particles_realloc(void *_ctx, int p, int new_n_particles)
{
  struct psc_mparticles *mprts = _ctx;
  struct psc_particles *prts = psc_mparticles_get_patch(mprts, p);
  struct psc_particles_cuda *cuda = psc_particles_cuda(prts);

  cuda->bnd_prts = realloc(cuda->bnd_prts, new_n_particles * sizeof(*cuda->bnd_prts));
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
// psc_bnd_particles_sub_setup

static void
psc_bnd_particles_sub_setup(struct psc_bnd_particles *bnd)
{
  bnd->ddcp = ddc_particles_create(bnd->psc->mrc_domain, sizeof(particle_t),
				   sizeof(particle_real_t),
				   MPI_PARTICLES_REAL,
				   ddcp_particles_realloc,
				   ddcp_particles_get_addr);
}

// ----------------------------------------------------------------------
// psc_bnd_particles_sub_unsetup

static void
psc_bnd_particles_sub_unsetup(struct psc_bnd_particles *bnd)
{
  ddc_particles_destroy(bnd->ddcp);
}

// ----------------------------------------------------------------------
// xchg_append helper

static void
xchg_append(struct psc_particles *prts, struct ddcp_patch *ddcp_patch, particle_t *prt)
{
  struct psc_particles_cuda *cuda = psc_particles_cuda(prts);
  cuda->bnd_prts[ddcp_patch->head++] = *prt;
}

static inline particle_t *
xchg_get_one(struct psc_particles *prts, int n)
{
  struct psc_particles_cuda *cuda = psc_particles_cuda(prts);
  return &cuda->bnd_prts[n];
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

static inline int
get_n_send(struct psc_particles *prts)
{
  struct psc_particles_cuda *cuda = psc_particles_cuda(prts);
  return cuda->bnd_n_send;
}

static inline int
get_head(struct psc_particles *prts)
{
  return 0;
}

#include "../psc_bnd_particles/psc_bnd_particles_exchange_particles_pre.c"

// ----------------------------------------------------------------------
// mprts_exchange_particles_pre

static void
mprts_exchange_particles_pre(struct psc_bnd_particles *bnd, struct psc_mparticles *mprts)
{
  for (int p = 0; p < mprts->nr_patches; p++) {
    struct psc_particles *prts = psc_mparticles_get_patch(mprts, p);
    exchange_particles_pre(bnd, prts);
  }
}

// ----------------------------------------------------------------------
// mprts_convert_to_cuda

static void
mprts_convert_to_cuda(struct psc_bnd_particles *bnd, struct psc_mparticles *mprts)
{
  struct psc_mparticles_cuda *mprts_cuda = psc_mparticles_cuda(mprts);
  struct cuda_mparticles *cmprts = mprts_cuda->cmprts;

  unsigned int nr_recv = 0;
  for (int p = 0; p < mprts->nr_patches; p++) {
    struct psc_particles *prts = psc_mparticles_get_patch(mprts, p);
    struct ddc_particles *ddcp = bnd->ddcp;
    struct ddcp_patch *patch = &ddcp->patches[prts->p];
    nr_recv += patch->head;
  }

  mprts_cuda->h_bnd_xi4  = malloc(nr_recv * sizeof(*mprts_cuda->h_bnd_xi4));
  mprts_cuda->h_bnd_pxi4 = malloc(nr_recv * sizeof(*mprts_cuda->h_bnd_pxi4));
  mprts_cuda->h_bnd_idx  = malloc(nr_recv * sizeof(*mprts_cuda->h_bnd_idx));
  mprts_cuda->h_bnd_off  = malloc(nr_recv * sizeof(*mprts_cuda->h_bnd_off));

  cuda_mparticles_zero_h_bnd_cnt(mprts);

  unsigned int off = 0;
  for (int p = 0; p < mprts->nr_patches; p++) {
    struct psc_particles *prts = psc_mparticles_get_patch(mprts, p);
    struct psc_particles_cuda *cuda = psc_particles_cuda(prts);
    struct ddc_particles *ddcp = bnd->ddcp;
    struct ddcp_patch *patch = &ddcp->patches[prts->p];
    int n_recv = patch->head;
    cuda->bnd_n_recv = n_recv;
    
    float4 *h_bnd_xi4 = mprts_cuda->h_bnd_xi4 + off;
    float4 *h_bnd_pxi4 = mprts_cuda->h_bnd_pxi4 + off;
    unsigned int *h_bnd_idx = mprts_cuda->h_bnd_idx + off;
    unsigned int *h_bnd_off = mprts_cuda->h_bnd_off + off;
    for (int n = 0; n < n_recv; n++) {
      particle_single_t *prt = &cuda->bnd_prts[n];
      h_bnd_xi4[n].x  = prt->xi;
      h_bnd_xi4[n].y  = prt->yi;
      h_bnd_xi4[n].z  = prt->zi;
      h_bnd_xi4[n].w  = cuda_int_as_float(prt->kind);
      h_bnd_pxi4[n].x = prt->pxi;
      h_bnd_pxi4[n].y = prt->pyi;
      h_bnd_pxi4[n].z = prt->pzi;
      h_bnd_pxi4[n].w = prt->qni_wni;

      int b_pos[3];
      for (int d = 0; d < 3; d++) {
	float *xi = &h_bnd_xi4[n].x;
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
      unsigned int b = (b_pos[2] * cuda->b_mx[1] + b_pos[1]) * cuda->b_mx[0] + b_pos[0];
      assert(b < cmprts->n_blocks_per_patch);
      b += p * cmprts->n_blocks_per_patch;
      h_bnd_idx[n] = b;
      h_bnd_off[n] = mprts_cuda->h_bnd_cnt[b]++;
    }
    free(cuda->bnd_prts);
    off += n_recv;
  }
}

// ----------------------------------------------------------------------
// psc_bnd_particles_sub_exchange_mprts_prep

static void
psc_bnd_particles_sub_exchange_mprts_prep(struct psc_bnd_particles *bnd,
				struct psc_mparticles *mprts)
{
  static int pr_A, pr_B, pr_D, pr_E, pr_F, pr_B0, pr_B1;
  if (!pr_A) {
    pr_A = prof_register("xchg_bidx", 1., 0, 0);
    pr_B0= prof_register("xchg_reduce", 1., 0, 0);
    pr_B1= prof_register("xchg_n_send", 1., 0, 0);
    pr_B = prof_register("xchg_scan_send", 1., 0, 0);
    pr_D = prof_register("xchg_from_dev", 1., 0, 0);
    pr_E = prof_register("xchg_cvt_from", 1., 0, 0);
    pr_F = prof_register("xchg_pre", 1., 0, 0);
  }

  //prof_start(pr_A);
  //cuda_mprts_find_block_keys(mprts);
  //prof_stop(pr_A);
  
  prof_start(pr_B0);
  cuda_mprts_spine_reduce(mprts);
  prof_stop(pr_B0);

  prof_start(pr_B1);
  cuda_mprts_find_n_send(mprts);
  prof_stop(pr_B1);

  prof_start(pr_B);
  cuda_mprts_scan_send_buf_total(mprts);
  prof_stop(pr_B);

  prof_start(pr_D);
  cuda_mprts_copy_from_dev(mprts);
  prof_stop(pr_D);
  
  prof_start(pr_E);
  cuda_mprts_convert_from_cuda(mprts);
  prof_stop(pr_E);
  
  prof_start(pr_F);
  mprts_exchange_particles_pre(bnd, mprts);
  prof_stop(pr_F);
}

// ----------------------------------------------------------------------
// psc_bnd_particles_sub_exchange_mprts_post

static void
psc_bnd_particles_sub_exchange_mprts_post(struct psc_bnd_particles *bnd,
					  struct psc_mparticles *mprts)
{
  struct psc_mparticles_cuda *mprts_cuda = psc_mparticles_cuda(mprts);

  static int pr_A, pr_B, pr_C, pr_D, pr_E, pr_D1;
  if (!pr_A) {
    pr_A = prof_register("xchg_cvt_to", 1., 0, 0);
    pr_B = prof_register("xchg_to_dev", 1., 0, 0);
    pr_C = prof_register("xchg_bidx", 1., 0, 0);
    pr_D = prof_register("xchg_sort", 1., 0, 0);
    pr_D1= prof_register("xchg_upd_off", 1., 0, 0);
    pr_E = prof_register("xchg_reorder", 1., 0, 0);
  }

  prof_start(pr_A);
  mprts_convert_to_cuda(bnd, mprts);
  prof_stop(pr_A);

  prof_start(pr_B);
  cuda_mprts_copy_to_dev(mprts);
  prof_stop(pr_B);

  prof_start(pr_C);
  cuda_mprts_find_block_indices_3(mprts);
  prof_stop(pr_C);
  
  prof_start(pr_D);
  cuda_mprts_sort(mprts);
  prof_stop(pr_D);

  prof_start(pr_D1);
  cuda_mprts_update_offsets(mprts);
  prof_stop(pr_D1);
  
#if 0
  prof_start(pr_E);
  cuda_mprts_reorder(mprts);
  mprts_cuda->need_reorder = false;
  prof_stop(pr_E);
  //  cuda_mprts_check_ordered_total(mprts);
#else
  mprts_cuda->need_reorder = true;
#endif
  
}

// ----------------------------------------------------------------------
// psc_bnd_particles_sub_exchange_particles_serial_periodic
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
psc_bnd_particles_sub_exchange_particles_serial_periodic(struct psc_bnd_particles *psc_bnd_particles,
						mparticles_cuda_t *particles)
{
  assert(0);
#if 0
  static int pr_F, pr_G, pr_H;
  if (!pr_F) {
    pr_F = prof_register("xchg_bidx_ids", 1., 0, 0);
    pr_G = prof_register("xchg_sort_pairs", 1., 0, 0);
    pr_H = prof_register("xchg_reorder_off", 1., 0, 0);
  }

  cuda_exchange_particles(0, psc_mparticles_get_patch(particles, 0));

  // sort
  for (int p = 0; p < particles->nr_patches; p++) {
    struct psc_particles *prts = psc_mparticles_get_patch(particles, p);
    struct psc_particles_cuda *cuda = psc_particles_cuda(prts);

    prof_start(pr_F);
    cuda_find_block_indices(prts, cuda->h_dev->bidx);
    prof_stop(pr_F);

    prof_start(pr_G);
    sort_pairs_device_2(cuda->sort_ctx, cuda->h_dev->bidx,
			cuda->h_dev->alt_ids,
			prts->n_part,
			cuda->h_dev->offsets);
    prof_stop(pr_G);

    prof_start(pr_H);
    cuda_reorder(prts, cuda->h_dev->alt_ids);
    prof_stop(pr_H);
  }
#endif
}

// ----------------------------------------------------------------------
// psc_bnd_particles_sub_exchange_particles_general

static void
psc_bnd_particles_sub_exchange_particles_general(struct psc_bnd_particles *bnd,
				       mparticles_cuda_t *particles)
{
  struct ddc_particles *ddcp = bnd->ddcp;

  static int pr_A, pr_B, pr_C;
  if (!pr_A) {
    pr_A = prof_register("xchg_prep", 1., 0, 0);
    pr_B = prof_register("xchg_comm", 1., 0, 0);
    pr_C = prof_register("xchg_post", 1., 0, 0);
  }
  
  prof_start(pr_A);
  psc_bnd_particles_sub_exchange_mprts_prep(bnd, particles);
  prof_stop(pr_A);

  prof_start(pr_B);
  ddc_particles_comm(ddcp, particles);
  prof_stop(pr_B);

  prof_start(pr_C);
  psc_bnd_particles_sub_exchange_mprts_post(bnd, particles);
  prof_stop(pr_C);
}

// ----------------------------------------------------------------------
// psc_bnd_particles_sub_exchange_particles

static void
psc_bnd_particles_sub_exchange_particles(struct psc_bnd_particles *bnd,
			       mparticles_base_t *particles_base)
{
  int size;
  MPI_Comm_size(psc_bnd_particles_comm(bnd), &size);

  // This function only makes sense if it's called for particles already being of cuda
  // type. We could call _get_cuda(), but that wouldn't be happy if some particles were
  // not in the right patch in the first place.

  assert(strcmp(psc_mparticles_type(particles_base), "cuda") == 0);
  mparticles_cuda_t *particles = particles_base;

  if (size == 1 && ppsc->nr_patches == 1 && // FIXME !!!
      ppsc->domain.bnd_fld_lo[0] == BND_FLD_PERIODIC &&
      ppsc->domain.bnd_fld_lo[1] == BND_FLD_PERIODIC &&
      ppsc->domain.bnd_fld_lo[2] == BND_FLD_PERIODIC) {
    psc_bnd_particles_sub_exchange_particles_serial_periodic(bnd, particles);
  } else {
    psc_bnd_particles_sub_exchange_particles_general(bnd, particles);
  }
}

// ======================================================================
// psc_bnd_particles: subclass "cuda"

struct psc_bnd_particles_ops psc_bnd_particles_cuda_ops = {
  .name                    = "cuda",
  .setup                   = psc_bnd_particles_sub_setup,
  .unsetup                 = psc_bnd_particles_sub_unsetup,
  .exchange_particles      = psc_bnd_particles_sub_exchange_particles,
  .exchange_mprts_prep     = psc_bnd_particles_sub_exchange_mprts_prep,
  .exchange_mprts_post     = psc_bnd_particles_sub_exchange_mprts_post,
};

