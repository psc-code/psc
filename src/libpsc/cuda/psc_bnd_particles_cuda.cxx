
#include "psc_bnd_particles_private.h"
#include "psc_particles_cuda.h"

#include "cuda_iface.h"

#include <mrc_profile.h>

#define DDCP_TYPE DDCP_TYPE_CUDA

#include "../psc_bnd_particles/ddc_particles_inc.c"

using mparticles_t = mparticles_cuda_t;

// ======================================================================
// psc_bnd_particles: subclass "cuda"

// ----------------------------------------------------------------------
// psc_bnd_particles_cuda_exchange_mprts_prep
  
template<>
void psc_bnd_particles_sub<mparticles_t>::exchange_mprts_prep(mparticles_t mprts)
{
  mprts->bnd_prep();
  
  for (int p = 0; p < ddcp->nr_patches; p++) {
    ddc_particles<mparticles_t>::patch *dpatch = &ddcp->patches[p];
    dpatch->m_buf = mprts->bnd_get_buffer(p);
    dpatch->m_begin = 0;
  }
}

// ----------------------------------------------------------------------
// psc_bnd_particles_cuda_exchange_mprts_post

template<>
void psc_bnd_particles_sub<mparticles_t>::exchange_mprts_post(mparticles_t mprts)
{
  mprts->bnd_post();
  
  for (int p = 0; p < ddcp->nr_patches; p++) {
    ddcp->patches[p].m_buf = NULL;
  }
};

// ----------------------------------------------------------------------
// psc_bnd_particles: subclass "cuda"

struct psc_bnd_particles_ops_cuda : psc_bnd_particles_ops {
  psc_bnd_particles_ops_cuda() {
    name                    = "cuda";
    size                    = sizeof(psc_bnd_particles_sub<mparticles_t>);
    setup                   = psc_bnd_particles_sub<mparticles_t>::setup;
    unsetup                 = psc_bnd_particles_sub<mparticles_t>::unsetup;
    exchange_particles      = psc_bnd_particles_sub<mparticles_t>::exchange_particles;
  }
} psc_bnd_particles_cuda_ops;

