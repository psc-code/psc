
#ifndef PSC_PARTICLE_C_H
#define PSC_PARTICLE_C_H

#include "psc.h"

typedef double particle_c_real_t;

#define MPI_PARTICLES_C_REAL MPI_DOUBLE

typedef struct {
  particle_c_real_t xi, yi, zi;
  particle_c_real_t pxi, pyi, pzi;
  particle_c_real_t qni;
  particle_c_real_t mni;
  particle_c_real_t wni;
} particle_c_t;

typedef struct {
  particle_c_t *particles;
  int n_part;
  int n_alloced;
} particles_c_t;

void particles_c_realloc(particles_c_t *pp, int new_n_part);

static inline particle_c_t *
particles_c_get_one(particles_c_t *pp, int n)
{
  return &pp->particles[n];
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

#endif
