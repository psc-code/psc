
#ifndef PSC_PARTICLE_C_H
#define PSC_PARTICLE_C_H

#include "psc.h"

#include <math.h>

typedef double particle_c_real_t;

#define MPI_PARTICLES_C_REAL MPI_DOUBLE

typedef struct psc_particle_c {
  particle_c_real_t xi, yi, zi;
  particle_c_real_t pxi, pyi, pzi;
  particle_c_real_t qni;
  particle_c_real_t mni;
  particle_c_real_t wni;
  long              kind; // 64 bits to match the other members, for bnd exchange
} particle_c_t;

typedef struct psc_particles_c {
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

static inline particle_c_real_t
particle_c_qni(particle_c_t *p)
{
  return p->qni;
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

static inline void
particle_c_get_relative_pos(particle_c_t *p, double xb[3],
			    particle_c_real_t xi[3])
{
  xi[0] = p->xi - xb[0];
  xi[1] = p->yi - xb[1];
  xi[2] = p->zi - xb[2];
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
