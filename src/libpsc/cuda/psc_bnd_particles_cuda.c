
#include "psc_bnd_particles_private.h"

#include "psc_bnd_cuda.h"
#include "psc_cuda.h"
#include "particles_cuda.h"
#include "psc_particles_as_cuda.h"
#include "cuda_mparticles.h"

#include <mrc_profile.h>

#define DDCP_TYPE DDCP_TYPE_CUDA

static void psc_bnd_particles_sub_exchange_mprts_prep_cuda(struct psc_bnd_particles *bnd,
							   struct psc_mparticles *mprts);
static void psc_bnd_particles_sub_exchange_mprts_post_cuda(struct psc_bnd_particles *bnd,
							   struct psc_mparticles *mprts);

#include "../psc_bnd_particles/ddc_particles_inc.c"

// ----------------------------------------------------------------------
// psc_bnd_particles_sub_exchange_mprts_prep_cuda

static void
psc_bnd_particles_sub_exchange_mprts_prep_cuda(struct psc_bnd_particles *bnd,
					       struct psc_mparticles *mprts)
{
  static int pr_A, pr_B, pr_D, pr_E, pr_B0, pr_B1;
  if (!pr_A) {
    pr_A = prof_register("xchg_bidx", 1., 0, 0);
    pr_B0= prof_register("xchg_reduce", 1., 0, 0);
    pr_B1= prof_register("xchg_n_send", 1., 0, 0);
    pr_B = prof_register("xchg_scan_send", 1., 0, 0);
    pr_D = prof_register("xchg_from_dev", 1., 0, 0);
    pr_E = prof_register("xchg_cvt_from", 1., 0, 0);
  }

  struct ddc_particles *ddcp = bnd->ddcp;
  struct psc_mparticles_cuda *mprts_cuda = psc_mparticles_cuda(mprts);
  struct cuda_mparticles *cmprts = mprts_cuda->cmprts;

  //prof_start(pr_A);
  //cuda_mprts_find_block_keys(mprts);
  //prof_stop(pr_A);
  
  prof_start(pr_B0);
  cuda_mparticles_spine_reduce(cmprts);
  prof_stop(pr_B0);

  prof_start(pr_B1);
  cuda_mparticles_find_n_send(cmprts);
  prof_stop(pr_B1);

  prof_start(pr_B);
  cuda_mparticles_scan_send_buf_total(cmprts);
  prof_stop(pr_B);

  prof_start(pr_D);
  cuda_mparticles_copy_from_dev(cmprts);
  prof_stop(pr_D);

  for (int p = 0; p < cmprts->n_patches; p++) {
    particle_buf_ctor(&cmprts->bnd.bpatch[p].buf);
    particle_buf_reserve(&cmprts->bnd.bpatch[p].buf, cmprts->bnd.bpatch[p].n_send);
    particle_buf_resize(&cmprts->bnd.bpatch[p].buf, cmprts->bnd.bpatch[p].n_send);

    struct ddcp_patch *dpatch = &ddcp->patches[p];
    dpatch->m_buf = &cmprts->bnd.bpatch[p].buf;
    dpatch->m_begin = 0;
  }

  // this will fill the buffers above
  prof_start(pr_E);
  cuda_mparticles_convert_from_cuda(cmprts);
  prof_stop(pr_E);
}

// ----------------------------------------------------------------------
// mprts_convert_to_cuda

static void
mprts_convert_to_cuda(struct psc_bnd_particles *bnd, struct psc_mparticles *mprts)
{
  struct ddc_particles *ddcp = bnd->ddcp;
  struct psc_mparticles_cuda *mprts_cuda = psc_mparticles_cuda(mprts);
  struct cuda_mparticles *cmprts = mprts_cuda->cmprts;

  unsigned int nr_recv = 0;
  for (int p = 0; p < mprts->nr_patches; p++) {
    struct ddc_particles *ddcp = bnd->ddcp;
    struct ddcp_patch *patch = &ddcp->patches[p];
    nr_recv += particle_buf_size(patch->m_buf);
  }

  cmprts->bnd.h_bnd_xi4  = malloc(nr_recv * sizeof(*cmprts->bnd.h_bnd_xi4));
  cmprts->bnd.h_bnd_pxi4 = malloc(nr_recv * sizeof(*cmprts->bnd.h_bnd_pxi4));
  cmprts->bnd.h_bnd_idx  = malloc(nr_recv * sizeof(*cmprts->bnd.h_bnd_idx));
  cmprts->bnd.h_bnd_off  = malloc(nr_recv * sizeof(*cmprts->bnd.h_bnd_off));

  memset(cmprts->bnd.h_bnd_cnt, 0,
	 cmprts->n_blocks * sizeof(*cmprts->bnd.h_bnd_cnt));

  unsigned int off = 0;
  for (int p = 0; p < mprts->nr_patches; p++) {
    struct ddcp_patch *patch = &ddcp->patches[p];
    int n_recv = particle_buf_size(patch->m_buf);
    cmprts->bnd.bpatch[p].n_recv = n_recv;
    
    float4 *h_bnd_xi4 = cmprts->bnd.h_bnd_xi4 + off;
    float4 *h_bnd_pxi4 = cmprts->bnd.h_bnd_pxi4 + off;
    unsigned int *h_bnd_idx = cmprts->bnd.h_bnd_idx + off;
    unsigned int *h_bnd_off = cmprts->bnd.h_bnd_off + off;
    for (int n = 0; n < n_recv; n++) {
      particle_t *prt = &cmprts->bnd.bpatch[p].buf.m_data[n];
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
	b_pos[d] = particle_real_fint(xi[d] * cmprts->b_dxi[d]);
	if (b_pos[d] < 0 || b_pos[d] >= cmprts->b_mx[d]) {
	  printf("!!! xi %g %g %g\n", xi[0], xi[1], xi[2]);
	  printf("!!! d %d xi4[n] %g biy %d // %d\n",
		 d, xi[d], b_pos[d], cmprts->b_mx[d]);
	  if (b_pos[d] < 0) {
	    xi[d] = 0.f;
	  } else {
	    xi[d] *= (1. - 1e-6);
	  }
	}
	b_pos[d] = particle_real_fint(xi[d] * cmprts->b_dxi[d]);
	assert(b_pos[d] >= 0 && b_pos[d] < cmprts->b_mx[d]);
      }
      unsigned int b = (b_pos[2] * cmprts->b_mx[1] + b_pos[1]) * cmprts->b_mx[0] + b_pos[0];
      assert(b < cmprts->n_blocks_per_patch);
      b += p * cmprts->n_blocks_per_patch;
      h_bnd_idx[n] = b;
      h_bnd_off[n] = cmprts->bnd.h_bnd_cnt[b]++;
    }
    off += n_recv;
  }
}

// ----------------------------------------------------------------------
// psc_bnd_particles_sub_exchange_mprts_post_cuda

static void
psc_bnd_particles_sub_exchange_mprts_post_cuda(struct psc_bnd_particles *bnd,
					       struct psc_mparticles *mprts)
{
  struct psc_mparticles_cuda *mprts_cuda = psc_mparticles_cuda(mprts);
  struct cuda_mparticles *cmprts = mprts_cuda->cmprts;

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
  unsigned int n_prts_by_patch[mprts->nr_patches];
  cuda_mparticles_get_size_all(cmprts, n_prts_by_patch);
  cuda_mprts_sort(mprts, (int *) n_prts_by_patch); // FIXME cast
  // FIXME, is this necessary, or doesn't update_offsets() do this, too?
  cuda_mparticles_resize_all(cmprts, n_prts_by_patch);
  prof_stop(pr_D);

  prof_start(pr_D1);
  cuda_mprts_update_offsets(mprts);
  prof_stop(pr_D1);
  
  prof_start(pr_E);
#if 0
  cuda_mparticles_reorder(cmprts);
  //  cuda_mprts_check_ordered_total(mprts);
#else
  cmprts->need_reorder = true;
#endif
  prof_stop(pr_E);

  struct ddc_particles *ddcp = bnd->ddcp;
  for (int p = 0; p < ddcp->nr_patches; p++) {
    struct ddcp_patch *dpatch = &ddcp->patches[p];
    particle_buf_dtor(&cmprts->bnd.bpatch[p].buf); // FIXME if we use it temporarily, it doesn't need to be in cmprts
    dpatch->m_buf = NULL;
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

