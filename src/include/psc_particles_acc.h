
#ifndef PSC_PARTICLE_ACC_H
#define PSC_PARTICLE_ACC_H

#include "psc_particles_private.h"

#include "cuda_wrap.h" // FIXME
#include <math.h>

typedef float particle_acc_real_t;


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


