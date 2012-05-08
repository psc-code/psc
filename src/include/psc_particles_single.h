
#ifndef PSC_PARTICLE_SINGLE_H
#define PSC_PARTICLE_SINGLE_H

#include "psc.h"

typedef float particle_single_real_t;

#define MPI_PARTICLES_SINGLE_REAL MPI_FLOAT

typedef struct {
  particle_single_real_t xi, yi, zi;
  particle_single_real_t pxi, pyi, pzi;
  particle_single_real_t qni_div_mni;
  particle_single_real_t qni_wni;
} particle_single_t;

typedef struct {
  particle_single_t *particles;
  int n_part;
  int n_alloced;
} particles_single_t;

void particles_single_realloc(particles_single_t *pp, int new_n_part);

static inline particle_single_t *
particles_single_get_one(particles_single_t *pp, int n)
{
  return &pp->particles[n];
}

static inline particle_single_real_t
particle_single_qni_div_mni(particle_single_t *p)
{
  return p->qni_div_mni;
}

static inline particle_single_real_t
particle_single_qni_wni(particle_single_t *p)
{
  return p->qni_wni;
}

static inline particle_single_real_t
particle_single_qni(particle_single_t *p)
{
  if (p->qni_wni > 0.f) {
    return 1.f;
  } else if (p->qni_wni < 0.f) {
    return -1.f;
  } else {
    assert(0);
  }
}

static inline particle_single_real_t
particle_single_wni(particle_single_t *p)
{
  return p->qni_wni / particle_single_qni(p);
}

static inline void
particle_single_get_relative_pos(particle_single_t *p, double xb[3],
				 particle_single_real_t xi[3])
{
  xi[0] = p->xi;
  xi[1] = p->yi;
  xi[2] = p->zi;
}

static inline int
particle_single_real_nint(particle_single_real_t x)
{
  return (int)(x + 10.5f) - 10;
}

static inline int
particle_single_real_fint(particle_single_real_t x)
{
  return (int)(x + 10.f) - 10;
}

#endif
