
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
} particles_c_t;

typedef struct {
  particles_c_t *p;
} mparticles_c_t;

void particles_c_alloc(particles_c_t *pp, int n_part);
void particles_c_realloc(particles_c_t *pp, int new_n_part);
void particles_c_free(particles_c_t *pp);
void particles_c_get(mparticles_c_t *particles, void *particles_base);
void particles_c_put(mparticles_c_t *particles, void *particles_base);

static inline particle_c_t *
particles_c_get_one(particles_c_t *pp, int n)
{
  return &pp->particles[n];
}

#endif
