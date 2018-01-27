
#include "psc_bnd_particles_private.h"
#include "psc_particles_cuda.h"

#include "cuda_iface.h"

#include <mrc_profile.h>

#define DDCP_TYPE DDCP_TYPE_CUDA

#include "../psc_bnd_particles/ddc_particles_inc.c"

using mparticles_t = mparticles_cuda_t;

// ======================================================================
// psc_bnd_particles: subclass "cuda"

template<typename MP>
struct bnd_particles_policy_cuda
{
  using mparticles_t = MP;
  using ddcp_t = ddc_particles<mparticles_t>;
  using ddcp_patch = typename ddcp_t::patch;
  
  // ----------------------------------------------------------------------
  // exchange_mprts_prep
  
  void exchange_mprts_prep(ddcp_t* ddcp, mparticles_t mprts)
  {
    mprts->bnd_prep();
    
    for (int p = 0; p < ddcp->nr_patches; p++) {
      ddcp_patch *dpatch = &ddcp->patches[p];
      dpatch->m_buf = mprts->bnd_get_buffer(p);
      dpatch->m_begin = 0;
    }
  }

  // ----------------------------------------------------------------------
  // exchange_mprts_post
  
  void exchange_mprts_post(ddcp_t* ddcp, mparticles_t mprts)
  {
    mprts->bnd_post();
    
    for (int p = 0; p < ddcp->nr_patches; p++) {
      ddcp->patches[p].m_buf = NULL;
    }
  }
};

// ----------------------------------------------------------------------
// psc_bnd_particles: subclass "cuda"

struct psc_bnd_particles_ops_cuda : psc_bnd_particles_ops {
  using sub_t = psc_bnd_particles_sub<mparticles_t,
				      bnd_particles_policy_cuda<mparticles_t>>;
  psc_bnd_particles_ops_cuda() {
    name                    = "cuda";
    size                    = sizeof(sub_t);
    setup                   = sub_t::setup;
    unsetup                 = sub_t::unsetup;
    exchange_particles      = sub_t::exchange_particles;
  }
} psc_bnd_particles_cuda_ops;

