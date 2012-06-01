
#ifndef PSC_PARTICLE_CUDA_H
#define PSC_PARTICLE_CUDA_H

#include "psc_particles_private.h"
#include "cuda_wrap.h"

#include "psc_particles_single.h"

typedef float particle_cuda_real_t;

#define MPI_PARTICLES_CUDA_REAL MPI_FLOAT

typedef struct {
  float4 *xi4;    // xi , yi , zi , qni_div_mni
  float4 *pxi4;   // pxi, pyi, pzi, qni_wni (if qni==0, then qni_wni = wni)
  float4 *alt_xi4;
  float4 *alt_pxi4;
  int *offsets;   // particles per block are
                  // are at indices offsets[block] .. offsets[block+1]-1
  int *c_offsets; // particles per cell offsets
  void *sort_ctx; // for sorting / particle xchg
  unsigned int *bidx;      // for particle xchg
  unsigned int *ids;       // for particle xchg
  unsigned int *alt_bidx;  // for particle xchg
  unsigned int *alt_ids;   // for particle xchg
  unsigned int *sums;      // for particle xchg
} particles_cuda_dev_t;

struct psc_particles_cuda {
  particles_cuda_dev_t d_part; // all particles, on device
  int nr_blocks;               // number of blocks
  int b_mx[3];                 // number of blocks by direction
  int n_alloced;
  int blocksize[3];            // dimensions of sub blocks in a patch
  particle_cuda_real_t b_dxi[3];
  struct cell_map map;         // maps 3d block pos to 1d block index

  // for bnd exchange
  particle_single_t *bnd_prts;
  float4 *bnd_xi4;
  float4 *bnd_pxi4;
  int bnd_n_part;
  int bnd_n_send;
  int bnd_n_part_save;
  unsigned int *bnd_cnt;
  unsigned int *bnd_idx;
  unsigned int *bnd_off;
};

#define psc_particles_cuda(prts) mrc_to_subobj(prts, struct psc_particles_cuda)

EXTERN_C void particles_cuda_get(struct psc_particles *pp);
EXTERN_C void particles_cuda_put(struct psc_particles *pp);

static inline int
particle_cuda_real_nint(particle_cuda_real_t x)
{
  return (int)(x + 10.5f) - 10;
}

static inline int
particle_cuda_real_fint(particle_cuda_real_t x)
{
  return (int)(x + 10.f) - 10;
}

#endif
