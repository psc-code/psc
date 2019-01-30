
#ifndef PSC_PARTICLE_SSE2_H
#define PSC_PARTICLE_SSE2_H

#include "psc.h"
#include "simd_sse2.h"

typedef sse2_real particle_sse2_real_t;

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

void particles_sse2_alloc(particles_sse2_t *pp, int n_part);
void particles_sse2_realloc(particles_sse2_t *pp, int new_n_part);
void particles_sse2_free(particles_sse2_t *pp);

static inline particle_sse2_t *
particles_sse2_get_one(particles_sse2_t *pp, int n)
{
  return &pp->particles[n];
}

#endif
