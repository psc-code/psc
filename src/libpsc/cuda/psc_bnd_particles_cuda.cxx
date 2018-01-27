
#include "psc_bnd_particles_private.h"

#include "psc_particles_as_cuda.h"
#include "cuda_iface.h"

#include <mrc_profile.h>

#define DDCP_TYPE DDCP_TYPE_CUDA

#include "../psc_bnd_particles/ddc_particles_inc.c"

// ======================================================================
// psc_bnd_particles: subclass "cuda"

template<>
struct mparticles_ddcp<mparticles_cuda_t>
{
  // ----------------------------------------------------------------------
  // psc_bnd_particles_cuda_exchange_mprts_prep
  
  static void exchange_mprts_prep(struct psc_bnd_particles *bnd,
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

  static void exchange_mprts_post(struct psc_bnd_particles *bnd,
				  mparticles_t mprts)
  {
    psc_bnd_particles_sub<mparticles_t>* sub = psc_bnd_particles_sub(bnd);

    mprts->bnd_post();

    ddc_particles<mparticles_t>* ddcp = static_cast<ddc_particles<mparticles_t>*>(sub->ddcp);
    for (int p = 0; p < ddcp->nr_patches; p++) {
      ddcp->patches[p].m_buf = NULL;
    }
  }

  // ----------------------------------------------------------------------
  // psc_bnd_particles_cuda_exchange_particles

  static void exchange_particles(struct psc_bnd_particles *bnd,
				 struct psc_mparticles *mprts_base)
  {
    // This function only makes sense if it's called for particles already being of cuda
    // type. If particles aren't in the right patches, the conversion in get_as would fail...

    assert(strcmp(psc_mparticles_type(mprts_base), "cuda") == 0);
    mparticles_t mprts = mprts_base->get_as<mparticles_t>();
  
    psc_bnd_particles_sub_exchange_particles_general<mparticles_cuda_t>(bnd, mprts);

    mprts.put_as(mprts_base);
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
    exchange_particles      = mparticles_ddcp<mparticles_cuda_t>::exchange_particles;
  }
} psc_bnd_particles_cuda_ops;

