
#ifndef PSC_PARTICLES_CBE_H
#define PSC_PARTICLES_CBE_H

#include "psc.h"
#include "simd_cbe.h"

typedef cbe_real particle_cbe_real_t;

#define MPI_PARTICLES_CBE_REAL MPI_CBE_REAL

#if CBE_DOUBLE

typedef struct {
  particle_cbe_real_t xi, yi, zi;
  particle_cbe_real_t pxi, pyi, pzi;
  particle_cbe_real_t qni, mni, wni;
  particle_cbe_real_t cni;
  // cni is just padding now. If I can managed to cut
  // one more element (possibly using the q/m ratio
  // modification) we can get rid of the padding all together. 
} particle_cbe_t;

#else

typedef struct {
  particle_cbe_real_t xi, yi, zi;
  particle_cbe_real_t pxi, pyi, pzi;
  particle_cbe_real_t qni, mni, wni;
  particle_cbe_real_t cni;
  particle_cbe_real_t pad1, pad2;
  // cni is just padding now. If I can managed to cut
  // one more element (possibly using the q/m ratio
  // modification) we can get rid of the padding all together. 
} particle_cbe_t;

#endif

typedef struct {
  particle_cbe_t *particles;
  int n_part;
} particles_cbe_t;

void particles_cbe_alloc(particles_cbe_t *pp, int n_part);
void particles_cbe_realloc(particles_cbe_t *pp, int new_n_part);
void particles_cbe_free(particles_cbe_t *pp);
void particles_cbe_get(particles_cbe_t *pp);
void particles_cbe_put(particles_cbe_t *pp);

static inline particle_cbe_t *
particles_cbe_get_one(particles_cbe_t *pp, int n)
{
  return &pp->particles[n];
}

#endif
