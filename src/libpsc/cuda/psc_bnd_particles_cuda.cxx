
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
  
template<>
void psc_bnd_particles_sub<mparticles_t>::exchange_mprts_prep(struct psc_bnd_particles *bnd,
							     mparticles_t mprts)
{
  psc_bnd_particles_sub<mparticles_t>* sub = psc_bnd_particles_sub(bnd);
  
  mprts->bnd_prep();
  
  ddc_particles<mparticles_t>* ddcp = static_cast<ddc_particles<mparticles_t>*>(sub->ddcp);
  for (int p = 0; p < ddcp->nr_patches; p++) {
    ddc_particles<mparticles_t>::patch *dpatch = &ddcp->patches[p];
    dpatch->m_buf = mprts->bnd_get_buffer(p);
    dpatch->m_begin = 0;
  }
}

// ----------------------------------------------------------------------
// psc_bnd_particles_cuda_exchange_mprts_post

template<>
void psc_bnd_particles_sub<mparticles_t>::exchange_mprts_post(struct psc_bnd_particles *bnd,
							      mparticles_t mprts)
{
  psc_bnd_particles_sub<mparticles_t>* sub = psc_bnd_particles_sub(bnd);
  
  mprts->bnd_post();
  
  ddc_particles<mparticles_t>* ddcp = static_cast<ddc_particles<mparticles_t>*>(sub->ddcp);
  for (int p = 0; p < ddcp->nr_patches; p++) {
    ddcp->patches[p].m_buf = NULL;
  }
};

// ----------------------------------------------------------------------
// psc_bnd_particles: subclass "cuda"

struct psc_bnd_particles_ops_cuda : psc_bnd_particles_ops {
  psc_bnd_particles_ops_cuda() {
    name                    = "cuda";
    size                    = sizeof(psc_bnd_particles_sub<mparticles_cuda_t>);
    setup                   = psc_bnd_particles_sub_setup;
    unsetup                 = psc_bnd_particles_sub_unsetup;
    exchange_particles      = psc_bnd_particles_sub_exchange_particles;
  }
} psc_bnd_particles_cuda_ops;

