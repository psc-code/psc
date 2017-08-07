
#ifndef PSC_PARTICLE_CUDA_H
#define PSC_PARTICLE_CUDA_H

#include "psc_particles_private.h"
#include "cuda_wrap.h"

#include "psc_particles_single.h"

typedef float particle_cuda_real_t;

#define MPI_PARTICLES_CUDA_REAL MPI_FLOAT

struct psc_particles_cuda {
  // for bnd exchange
  particle_single_t *bnd_prts;
  int bnd_n_recv;
  int bnd_n_send;
  void *sort_ctx; // for sorting / particle xchg
};

#define psc_particles_cuda(prts) mrc_to_subobj(prts, struct psc_particles_cuda)

struct psc_mparticles_cuda {
  struct cuda_mparticles *cmprts;

  unsigned int *d_alt_bidx;
  unsigned int *d_sums; // FIXME, too many arrays, consolidation would be good
  unsigned int nr_prts_send;
  unsigned int nr_prts_recv;
  unsigned int *d_bnd_spine_cnts;
  unsigned int *d_bnd_spine_sums;
  float4 *h_bnd_xi4, *h_bnd_pxi4;
  unsigned int *h_bnd_idx;
  unsigned int *h_bnd_off;
  unsigned int *h_bnd_cnt;

  int b_mx[3];                 // number of blocks by direction
  particle_cuda_real_t b_dxi[3];
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
