
#ifndef CUDA_MPARTICLES_H
#define CUDA_MPARTICLES_H

struct cuda_mparticles {
  float4 *d_xi4, *d_pxi4;         // current particle data
  float4 *d_alt_xi4, *d_alt_pxi4; // storage for out-of-place reordering of particle data
};
  
#endif
