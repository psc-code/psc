
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
} psc_particles_c_t;

static inline particle_c_t *
psc_particles_c_get_one(psc_particles_c_t *pp, int n)
{
  return &pp->particles[n];
}

#endif
