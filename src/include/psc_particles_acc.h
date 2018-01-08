
#ifndef PSC_PARTICLE_ACC_H
#define PSC_PARTICLE_ACC_H

#include "psc_particles_private.h"

#include "cuda_wrap.h" // FIXME
#include <math.h>

typedef float particle_acc_real_t;


#define MPI_PARTICLES_ACC_REAL MPI_FLOAT

typedef struct psc_particle_acc {
  particle_acc_real_t xi[3];
  particle_acc_real_t kind_as_float;
  particle_acc_real_t pxi[3];
  particle_acc_real_t qni_wni;
} particle_acc_t;

struct psc_particles_acc {
  float4 *xi4, *pxi4;
  float4 *xi4_alt, *pxi4_alt;
  unsigned int *b_idx;
  unsigned int *b_ids;
  unsigned int *b_cnt;
  unsigned int *b_off;

  particle_acc_real_t dxi[3];
  int b_mx[3];
  int nr_blocks;
};

#define psc_particles_acc(prts) mrc_to_subobj(prts, struct psc_particles_acc)

struct psc_mparticles_acc {
  int n_part_total; // total in all patches
  int n_alloced_total;

  particle_acc_real_t dxi[3];
  int b_mx[3];
  int bs[3];
  int nr_blocks;
  int nr_blocks_total;

  float4 *xi4, *pxi4;
  float4 *xi4_alt, *pxi4_alt;
  unsigned int *b_off;
};

#define psc_mparticles_acc(prts) mrc_to_subobj(prts, struct psc_mparticles_acc)

#define particle_acc_wni(p) ({				\
      particle_acc_real_t rv;					\
      int kind = particle_acc_kind(p);			\
      rv = p->qni_wni / ppsc->kinds[kind].q;			\
      rv;							\
    })

static inline int
particle_acc_kind(particle_acc_t *prt)
{
  return cuda_float_as_int(prt->kind_as_float);
}

#define particle_acc_qni_wni(prt) ((prt)->qni_wni)

#define particle_acc_px(prt) ((prt)->pxi[0])

static inline int
particle_acc_real_fint(particle_acc_real_t x)
{
#ifdef __CUDACC__
  return __float2int_rd(x);
#else
  return floorf(x);
#endif
}

static inline particle_acc_real_t
particle_acc_real_sqrt(particle_acc_real_t x)
{
  return sqrtf(x);
}

static inline particle_acc_real_t
particle_acc_real_abs(particle_acc_real_t x)
{
  return fabsf(x);
}

#define PARTICLE_ACC_LOAD_POS(prt, d_xi4, n) do {			\
    float4 xi4 = d_xi4[n];						\
    (prt).xi[0] = xi4.x;						\
    (prt).xi[1] = xi4.y;						\
    (prt).xi[2] = xi4.z;						\
    (prt).kind_as_float = xi4.w;					\
  } while (0)

#define PARTICLE_ACC_LOAD_MOM(prt, d_pxi4, n) do {			\
    float4 pxi4 = d_pxi4[n];						\
    (prt).pxi[0]  = pxi4.x;						\
    (prt).pxi[1]  = pxi4.y;						\
    (prt).pxi[2]  = pxi4.z;						\
    (prt).qni_wni = pxi4.w;						\
  } while (0)

#define PARTICLE_ACC_STORE_POS(pp, d_xi4, n) do {			\
    float4 xi4 = { (pp).xi[0], (pp).xi[1], (pp).xi[2], (pp).kind_as_float }; \
    d_xi4[n] = xi4;							\
} while (0)

#define PARTICLE_ACC_STORE_MOM(pp, d_pxi4, n) do {			\
    float4 pxi4 = { (pp).pxi[0], (pp).pxi[1], (pp).pxi[2], (pp).qni_wni }; \
    d_pxi4[n] = pxi4;							\
} while (0)


#endif


