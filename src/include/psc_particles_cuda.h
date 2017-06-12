
#ifndef PSC_PARTICLE_CUDA_H
#define PSC_PARTICLE_CUDA_H

#include "psc_particles_private.h"
#include "cuda_wrap.h"

#include "psc_particles_single.h"

typedef float particle_cuda_real_t;

#define MPI_PARTICLES_CUDA_REAL MPI_FLOAT

typedef struct {
  int n_part;     // # of particles in this patch
} particles_cuda_dev_t;

struct psc_particles_cuda {
  particles_cuda_dev_t *h_dev; // info that we also keep on the device, but this is on the host
  particles_cuda_dev_t *d_dev; // this one actually lives in device mem
  int nr_blocks;               // number of blocks
  int b_mx[3];                 // number of blocks by direction
  int blocksize[3];            // dimensions of sub blocks in a patch
  particle_cuda_real_t b_dxi[3];
  struct cell_map map;         // maps 3d block pos to 1d block index

  // for bnd exchange
  particle_single_t *bnd_prts;
  int bnd_n_recv;
  int bnd_n_send;
  void *sort_ctx; // for sorting / particle xchg
  struct psc_mparticles *mprts; // parent containing this patch of particles
};

#define psc_particles_cuda(prts) mrc_to_subobj(prts, struct psc_particles_cuda)

struct psc_mparticles_cuda {
  struct cuda_mparticles *cmprts;
  particles_cuda_dev_t *h_dev; // particles_cuda_dev_t array for all patches
  particles_cuda_dev_t *d_dev; // same, living in device memory
  unsigned int nr_prts;     // total # of particles in all patches
  unsigned int nr_alloced;  // arrays are alloced for this # of particles

  float4 *d_alt_xi4, *d_alt_pxi4;
  unsigned int *d_bidx;
  unsigned int *d_alt_bidx;
  unsigned int *d_ids;
  unsigned int *d_sums; // FIXME, too many arrays, consolidation would be good
  unsigned int nr_prts_send;
  unsigned int nr_prts_recv;
  unsigned int *d_off; // particles per block
                       // are at indices offsets[block] .. offsets[block+1]-1
                       // indices numbered for the total mprts array 
  unsigned int *d_bnd_spine_cnts;
  unsigned int *d_bnd_spine_sums;
  float4 *h_bnd_xi4, *h_bnd_pxi4;
  unsigned int *h_bnd_idx;
  unsigned int *h_bnd_off;
  unsigned int *h_bnd_cnt;
  int nr_blocks;                 // number of blocks per patch
  int nr_total_blocks;           // number of blocks in all patches in mprts
  int b_mx[3];                   // number of blocks by direction
  int blocksize[3];              // dimensions of sub blocks in a patch
  particle_cuda_real_t b_dxi[3]; // 1. / (blocksize[d] * dx[d])
  bool need_reorder;
};

#define psc_mparticles_cuda(prts) mrc_to_subobj(prts, struct psc_mparticles_cuda)

#define CUDA_BND_S_NEW (9)
#define CUDA_BND_S_OOB (10)
#define CUDA_BND_STRIDE (10)

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
