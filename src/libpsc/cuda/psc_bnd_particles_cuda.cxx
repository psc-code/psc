
#include "psc_bnd_particles_private.h"

#include "psc_particles_as_cuda.h"
#include "cuda_iface.h"

#include <mrc_profile.h>

#define DDCP_TYPE DDCP_TYPE_CUDA

#include "../psc_bnd_particles/ddc_particles_inc.c"

// ======================================================================
// psc_bnd_particles: subclass "cuda"

// ----------------------------------------------------------------------
// psc_bnd_particles_cuda_exchange_mprts_prep

static void
psc_bnd_particles_cuda_exchange_mprts_prep(struct psc_bnd_particles *bnd,
					   struct psc_mparticles *_mprts)
{
  mparticles_cuda_t mprts(_mprts);
  
  mprts->bnd_prep();
  
  ddc_particles<mparticles_t>* ddcp = static_cast<ddc_particles<mparticles_t>*>(bnd->ddcp);
  for (int p = 0; p < ddcp->nr_patches; p++) {
    ddc_particles<mparticles_t>::patch *dpatch = &ddcp->patches[p];
    dpatch->m_buf = mprts->bnd_get_buffer(p);
    dpatch->m_begin = 0;
  }
}

// ----------------------------------------------------------------------
// psc_bnd_particles_cuda_exchange_mprts_post

static void
psc_bnd_particles_cuda_exchange_mprts_post(struct psc_bnd_particles *bnd,
					   struct psc_mparticles *_mprts)
{
  mparticles_t mprts(_mprts);

  mprts->bnd_post();

  ddc_particles<mparticles_t>* ddcp = static_cast<ddc_particles<mparticles_t>*>(bnd->ddcp);
  for (int p = 0; p < ddcp->nr_patches; p++) {
    ddcp->patches[p].m_buf = NULL;
  }
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
//
// or maybe (probably) this special case should just be killed, since we don't
// care about the single GPU case that much...

static void
psc_bnd_particles_sub_exchange_particles_serial_periodic(struct psc_bnd_particles *psc_bnd_particles,
						struct psc_mparticles *mprts)
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
    prof_start(pr_F);
    cuda_find_block_indices(prts, cuda->h_dev->bidx);
    prof_stop(pr_F);

    prof_start(pr_G);
    sort_pairs_device_2(cuda->sort_ctx, cuda->h_dev->bidx,
			cuda->h_dev->alt_ids,
			mprts_cuda->n_prts_by_patch[p],
			cuda->h_dev->offsets);
    prof_stop(pr_G);

    prof_start(pr_H);
    cuda_reorder(prts, cuda->h_dev->alt_ids);
    prof_stop(pr_H);
  }
#endif
}

// ----------------------------------------------------------------------
// psc_bnd_particles_cuda_exchange_particles

static void
psc_bnd_particles_cuda_exchange_particles(struct psc_bnd_particles *bnd,
					  struct psc_mparticles *mprts_base)
{
  // This function only makes sense if it's called for particles already being of cuda
  // type. If particles aren't in the right patches, the conversion in get_as would fail...

  assert(strcmp(psc_mparticles_type(mprts_base), "cuda") == 0);
  mparticles_t mprts = mprts_base->get_as<mparticles_t>();
  
  int size;
  MPI_Comm_size(psc_bnd_particles_comm(bnd), &size);

  if (size == 1 && ppsc->nr_patches == 1 && // FIXME !!!
      ppsc->domain.bnd_fld_lo[0] == BND_FLD_PERIODIC &&
      ppsc->domain.bnd_fld_lo[1] == BND_FLD_PERIODIC &&
      ppsc->domain.bnd_fld_lo[2] == BND_FLD_PERIODIC) {
    psc_bnd_particles_sub_exchange_particles_serial_periodic(bnd, mprts.mprts());
  } else {
    psc_bnd_particles_sub_exchange_particles_general(bnd, mprts.mprts());
  }

  mprts.put_as(mprts_base);
}

// ----------------------------------------------------------------------
// psc_bnd_particles: subclass "cuda"

struct psc_bnd_particles_ops_cuda : psc_bnd_particles_ops {
  psc_bnd_particles_ops_cuda() {
    name                    = "cuda";
    setup                   = psc_bnd_particles_sub_setup;
    unsetup                 = psc_bnd_particles_sub_unsetup;
    exchange_particles      = psc_bnd_particles_cuda_exchange_particles;
    exchange_mprts_prep     = psc_bnd_particles_cuda_exchange_mprts_prep;
    exchange_mprts_post     = psc_bnd_particles_cuda_exchange_mprts_post;
  }
} psc_bnd_particles_cuda_ops;

