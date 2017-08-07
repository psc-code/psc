
#ifndef PSC_PARTICLE_C_H
#define PSC_PARTICLE_C_H

#include "psc_particles_private.h"

#include <math.h>

typedef double particle_c_real_t;

#define MPI_PARTICLES_C_REAL MPI_DOUBLE

typedef struct psc_particle_c {
  particle_c_real_t xi, yi, zi;
  particle_c_real_t pxi, pyi, pzi;
  particle_c_real_t qni;
  particle_c_real_t mni;
  particle_c_real_t wni;
  long long         kind; // 64 bits to match the other members, for bnd exchange
} particle_c_t;

#define psc_particles_c(prts) mrc_to_subobj(prts, struct psc_particles_c)

struct psc_mparticles_c_patch {
  particle_c_t *prt_array;
};

struct psc_mparticles_c {
  struct psc_mparticles_c_patch *patch;
};

#define psc_mparticles_c(prts) mrc_to_subobj(prts, struct psc_mparticles_c)

static inline particle_c_t *
particles_c_get_one(struct psc_particles *prts, int n)
{
  struct psc_mparticles *mprts = prts->mprts;
  int p = prts->p;
  
  assert(psc_mparticles_ops(mprts) == &psc_mparticles_c_ops);
  return &psc_mparticles_c(mprts)->patch[p].prt_array[n];
}

static inline particle_c_real_t
particle_c_qni_div_mni(particle_c_t *p)
{
  return p->qni / p->mni;
}

static inline particle_c_real_t
particle_c_qni_wni(particle_c_t *p)
{
  return p->qni * p->wni;
}

static inline particle_c_real_t
particle_c_qni(particle_c_t *p)
{
  return p->qni;
}

static inline particle_c_real_t
particle_c_mni(particle_c_t *p)
{
  return p->mni;
}

static inline particle_c_real_t
particle_c_wni(particle_c_t *p)
{
  return p->wni;
}

static inline int
particle_c_kind(particle_c_t *p)
{
  return p->kind;
}

#define particle_c_x(prt) ((prt)->xi)
#define particle_c_px(prt) ((prt)->pxi)

static inline void
particle_c_get_relative_pos(particle_c_t *p, double xb[3],
			    particle_c_real_t xi[3])
{
  xi[0] = p->xi;
  xi[1] = p->yi;
  xi[2] = p->zi;
}

static inline int
particle_c_real_nint(particle_c_real_t x)
{
  return (int)(x + 10.5f) - 10;
}

static inline int
particle_c_real_fint(particle_c_real_t x)
{
  return (int)(x + 10.f) - 10;
}

static inline particle_c_real_t
particle_c_real_sqrt(particle_c_real_t x)
{
  return sqrt(x);
}

static inline particle_c_real_t
particle_c_real_abs(particle_c_real_t x)
{
  return fabs(x);
}

#endif
