
#ifndef PSC_PARTICLE_CUDA2_H
#define PSC_PARTICLE_CUDA2_H

#include "psc_particles_private.h"

typedef float particle_cuda2_real_t;

#define MPI_PARTICLES_CUDA2_REAL MPI_FLOAT

typedef struct psc_particle_cuda2 {
  particle_cuda2_real_t xi, yi, zi;
  particle_cuda2_real_t qni_wni;
  particle_cuda2_real_t pxi, pyi, pzi;
  int kind;
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
      rv = p->qni_wni / ppsc->kinds[p->kind].q;			\
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


