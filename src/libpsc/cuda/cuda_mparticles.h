
#ifndef CUDA_MPARTICLES_H
#define CUDA_MPARTICLES_H

struct cuda_mparticles {
  // per particle
  float4 *d_xi4, *d_pxi4;         // current particle data
  float4 *d_alt_xi4, *d_alt_pxi4; // storage for out-of-place reordering of particle data
  unsigned int *d_bidx;           // block index (incl patch) per particle
  unsigned int *d_id;             // particle id for sorting

  int n_patches;                  // # of patches
  unsigned int n_prts;            // total # of particles across all patches
  unsigned int n_alloced;         // size of particle-related arrays as allocated
};

struct cuda_mparticles *cuda_mparticles_create(void);
void cuda_mparticles_destroy(struct cuda_mparticles *cmprts);
void cuda_mparticles_alloc(struct cuda_mparticles *cmprts, int n_alloced);
void cuda_mparticles_dealloc(struct cuda_mparticles *cmprts);
void cuda_mparticles_swap_alt(struct cuda_mparticles *cmprts);

#endif
