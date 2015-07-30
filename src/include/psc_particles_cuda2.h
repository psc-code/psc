
#ifndef PSC_PARTICLE_CUDA2_H
#define PSC_PARTICLE_CUDA2_H

#include "psc_particles_private.h"

#include "cuda_wrap.h"

typedef float particle_cuda2_real_t;

#define MPI_PARTICLES_CUDA2_REAL MPI_FLOAT

struct psc_particles_cuda2 {
  // on host
  float4 *h_xi4, *h_pxi4;
  float4 *h_xi4_alt, *h_pxi4_alt;
  unsigned int *b_idx;
  unsigned int *b_ids;
  unsigned int *b_cnt;
  unsigned int *b_off;

  // on device
  float4 *d_xi4, *d_pxi4;
  unsigned int *d_b_off;

  int n_alloced;
  particle_cuda2_real_t dxi[3];
  int b_mx[3];
  int nr_blocks;
};

#define psc_particles_cuda2(prts) mrc_to_subobj(prts, struct psc_particles_cuda2)

struct psc_mparticles_cuda2 {
};

#define psc_mparticles_cuda2(prts) mrc_to_subobj(prts, struct psc_mparticles_cuda2)

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

#define _LOAD_PARTICLE_POS(prt, d_xi4, n) do {				\
    float4 xi4 = d_xi4[n];						\
    (prt).xi4.x = xi4.x;						\
    (prt).xi4.y = xi4.y;						\
    (prt).xi4.z = xi4.z;						\
    (prt).xi4.w = xi4.w;						\
  } while (0)

#define _LOAD_PARTICLE_MOM(prt, d_pxi4, n) do {				\
    float4 pxi4 = d_pxi4[n];						\
    (prt).pxi4.x = pxi4.x;						\
    (prt).pxi4.y = pxi4.y;						\
    (prt).pxi4.z = pxi4.z;						\
    (prt).pxi4.w = pxi4.w;						\
  } while (0)

#define _STORE_PARTICLE_POS(prt, d_xi4, n) do {				\
    float4 xi4;								\
    xi4.x = (prt).xi4.x;						\
    xi4.y = (prt).xi4.y;						\
    xi4.z = (prt).xi4.z;						\
    xi4.w = (prt).xi4.w;						\
    d_xi4[n] = xi4;							\
  } while (0)

#define _STORE_PARTICLE_MOM(prt, d_pxi4, n) do {			\
    float4 pxi4;							\
    pxi4.x = (prt).pxi4.x;						\
    pxi4.y = (prt).pxi4.y;						\
    pxi4.z = (prt).pxi4.z;						\
    pxi4.w = (prt).pxi4.w;						\
    d_pxi4[n] = pxi4;							\
  } while (0)

typedef struct psc_particle_cuda2 {
  float4 xi4;
  float4 pxi4;
} particle_cuda2_t;

#endif


