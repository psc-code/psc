
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
  struct cuda_mparticles *cmprts = psc_mparticles_cuda(mprts)->cmprts;

  cuda_mparticles_bnd_prep(cmprts);
  
  struct ddc_particles *ddcp = bnd->ddcp;
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

  cuda_mparticles_bnd_post(cmprts);

  struct ddc_particles *ddcp = bnd->ddcp;
  for (int p = 0; p < cmprts->n_patches; p++) {
    ddcp->patches[p].m_buf = NULL;
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

