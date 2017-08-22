
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

  prof_start(pr_E);
  cuda_mparticles_convert_from_cuda(cmprts);
  prof_stop(pr_E);

  for (int p = 0; p < cmprts->n_patches; p++) {
    struct ddcp_patch *dpatch = &ddcp->patches[p];
    dpatch->m_buf = &cmprts->bnd.bpatch[p].buf;
    dpatch->m_begin = 0;
  }
}

// ----------------------------------------------------------------------
// psc_bnd_particles_sub_exchange_mprts_post_cuda

static void
psc_bnd_particles_sub_exchange_mprts_post_cuda(struct psc_bnd_particles *bnd,
					       struct psc_mparticles *mprts)
{
  struct cuda_mparticles *cmprts = psc_mparticles_cuda(mprts)->cmprts;

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
  cuda_mparticles_convert_to_cuda(cmprts);
  prof_stop(pr_A);

  prof_start(pr_B);
  cuda_mparticles_copy_to_dev(cmprts);
  prof_stop(pr_B);

  prof_start(pr_C);
  cuda_mparticles_find_block_indices_3(cmprts);
  prof_stop(pr_C);
  
  prof_start(pr_D);
  unsigned int n_prts_by_patch[cmprts->n_patches];
  cuda_mparticles_get_size_all(cmprts, n_prts_by_patch);
  cuda_mparticles_sort(cmprts, (int *) n_prts_by_patch); // FIXME cast
  // FIXME, is this necessary, or doesn't update_offsets() do this, too?
  cuda_mparticles_resize_all(cmprts, n_prts_by_patch);
  prof_stop(pr_D);

  prof_start(pr_D1);
  cuda_mparticles_update_offsets(cmprts);
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
    particle_buf_dtor(&cmprts->bnd.bpatch[p].buf); // FIXME if we use it temporarily, it doesn't need to be in cmprts
    struct ddcp_patch *dpatch = &ddcp->patches[p];
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

