
#ifndef PSC_PARTICLE_CUDA_H
#define PSC_PARTICLE_CUDA_H

#include "psc.h"
#include "cuda_wrap.h"

typedef float particle_cuda_real_t;

#define MPI_PARTICLES_CUDA_REAL MPI_FLOAT

typedef struct {
  float4 *xi4;    // xi , yi , zi , qni_div_mni
  float4 *pxi4;   // pxi, pyi, pzi, qni_wni (if qni==0, then qni_wni = wni)
  int *offsets;   // particles per block are
                  // are at indices offsets[block] .. offsets[block+1]-1
} particles_cuda_dev_t;

typedef struct {
  particles_cuda_dev_t h_part; // all particles, on host
  particles_cuda_dev_t d_part; // all particles, on device
  int nr_blocks;               // number of blocks
  int b_mx[3];                 // number of blocks by direction
} particles_cuda_t;

//void particles_cuda_alloc(particles_cuda_t *pp, int n_part);
//void particles_cuda_realloc(particles_cuda_t *pp, int new_n_part);
//void particles_cuda_free(particles_cuda_t *pp);
EXTERN_C void particles_cuda_get(particles_cuda_t *pp);
EXTERN_C void particles_cuda_put(particles_cuda_t *pp);

//static inline particle_cuda_t *
//particles_cuda_get_one(particles_cuda_t *pp, int n)
//{
//  return &pp->particles[n];
//}

#endif
