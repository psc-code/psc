
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
} psc_particles_sse2_t;

//void psc_particles_sse2_alloc(struct psc_particles_c *pp, int n_part);
//void psc_particles_sse2_realloc(struct psc_particles_c *pp, int new_n_part);
//void psc_particles_sse2_free(struct psc_particles_c *pp);
void psc_particles_sse2_get(psc_particles_sse2_t *pp);
void psc_particles_sse2_put(psc_particles_sse2_t *pp);

static inline particle_sse2_t *
psc_particles_sse2_get_one(psc_particles_sse2_t *pp, int n)
{
  return &pp->particles[n];
}

#endif
