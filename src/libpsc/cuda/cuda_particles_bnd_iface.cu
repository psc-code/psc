
#include "cuda_particles_bnd_iface.h"
#include "cuda_mparticles.h"

// ----------------------------------------------------------------------
// prep

void cuda_particles_bnd::prep(ddcp_t* ddcp, cuda_mparticles* cmprts)
{
  cmprts->bnd_prep();

  for (int p = 0; p < ddcp->nr_patches; p++) {
    ddcp_patch *dpatch = &ddcp->patches[p];
    dpatch->m_buf = cmprts->bnd_get_buffer(p);
    dpatch->m_begin = 0;
  }
}

// ----------------------------------------------------------------------
// post

void cuda_particles_bnd::post(ddcp_t* ddcp, cuda_mparticles* cmprts)
{
  cmprts->bnd_post();
    
  for (int p = 0; p < ddcp->nr_patches; p++) {
    ddcp->patches[p].m_buf = NULL;
  }
}


