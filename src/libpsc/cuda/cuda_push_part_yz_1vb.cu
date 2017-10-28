
#include "cuda_mparticles.h"
#include "particles_cuda.h"

void
cuda_mparticles_params_set(struct cuda_mparticles_params *mprts_prm,
			   struct cuda_mparticles *cmprts)
{
  mprts_prm->fnqs = cmprts->fnqs;

  for (int d = 0; d < 3; d++) {
    mprts_prm->b_mx[d] = cmprts->b_mx[d];
    mprts_prm->b_dxi[d] = cmprts->b_dxi[d];
    mprts_prm->dxi[d] = 1.f / cmprts->dx[d];
  }
}

