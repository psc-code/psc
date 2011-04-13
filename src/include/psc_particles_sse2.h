
#ifndef PSC_PARTICLE_SSE2_H
#define PSC_PARTICLE_SSE2_H

#include "psc.h"
#include "simd_sse2.h"

typedef sse2_real particle_sse2_real_t;

#define MPI_PARTICLES_SSE2_REAL MPI_SSE2_REAL

typedef struct {
  particle_sse2_real_t xi, yi, zi;
  particle_sse2_real_t pxi, pyi, pzi;
  particle_sse2_real_t qni;
  particle_sse2_real_t mni;
  particle_sse2_real_t wni;
} particle_sse2_t;

typedef struct {
  particle_sse2_t *particles;
  int n_part;
} particles_sse2_t;

typedef struct {
  particles_sse2_t *p;
} mparticles_sse2_t;

void particles_sse2_alloc(particles_sse2_t *pp, int n_part);
void particles_sse2_realloc(particles_sse2_t *pp, int new_n_part);
void particles_sse2_free(particles_sse2_t *pp);
void particles_sse2_get(particles_sse2_t *pp, void *particles_base);
void particles_sse2_put(particles_sse2_t *pp, void *particles_base);

static inline particle_sse2_t *
particles_sse2_get_one(particles_sse2_t *pp, int n)
{
  return &pp->particles[n];
}

#endif
