
#ifndef PSC_PARTICLE_CUDA2_H
#define PSC_PARTICLE_CUDA2_H

#include "psc_particles_private.h"

#include "cuda_wrap.h"

typedef float particle_cuda2_real_t;

#define MPI_PARTICLES_CUDA2_REAL MPI_FLOAT

typedef struct psc_particle_cuda2 {
  float4 xi4;
  float4 pxi4;
} particle_cuda2_t;

struct psc_particles_cuda2 {
  particle_cuda2_t *particles;
  particle_cuda2_t *particles_alt;
  int n_alloced;
  particle_cuda2_real_t dxi[3];
  int b_mx[3];
  int nr_blocks;
  unsigned int *b_idx;
  unsigned int *b_ids;
  unsigned int *b_cnt;
  unsigned int *b_off;
};

#define psc_particles_cuda2(prts) mrc_to_subobj(prts, struct psc_particles_cuda2)

static inline particle_cuda2_t *
particles_cuda2_get_one(struct psc_particles *prts, int n)
{
  assert(psc_particles_ops(prts) == &psc_particles_cuda2_ops);
  return &psc_particles_cuda2(prts)->particles[n];
}

#define particle_cuda2_wni(p) ({				\
      particle_cuda2_real_t rv;					\
      int kind = cuda_float_as_int(p->xi4.w);			\
      rv = p->pxi4.w / ppsc->kinds[kind].q;			\
      rv;							\
    })

static inline int
particle_cuda2_real_fint(particle_cuda2_real_t x)
{
  return floorf(x);
}

static inline particle_cuda2_real_t
particle_cuda2_real_sqrt(particle_cuda2_real_t x)
{
  return sqrtf(x);
}

static inline particle_cuda2_real_t
particle_cuda2_real_abs(particle_cuda2_real_t x)
{
  return fabsf(x);
}

#endif


