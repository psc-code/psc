
#ifndef CUDA_MPARTICLES_H
#define CUDA_MPARTICLES_H

struct cuda_mparticles {
  float4 *d_xi4, *d_pxi4;         // current particle data
  float4 *d_alt_xi4, *d_alt_pxi4; // storage for out-of-place reordering of particle data

  int n_patches;
};

struct cuda_mparticles *cuda_mparticles_create(void);
void cuda_mparticles_destroy(struct cuda_mparticles *cmprts);
void cuda_mparticles_alloc(struct cuda_mparticles *cmprts, int nr_alloced);
void cuda_mparticles_dealloc(struct cuda_mparticles *cmprts);
void cuda_mparticles_swap_alt(struct cuda_mparticles *cmprts);

#endif
