#ifndef PSC_PARTICLES_CBE_H
#define PSC_PARTICLES_CBE_H

#include "simd_cbe.h"

typedef cbe_real particle_cbe_real_t;

#if CBE_DOUBLE

typedef struct psc_particle_cbe {
  particle_cbe_real_t xi, yi, zi;
  particle_cbe_real_t pxi, pyi, pzi;
  particle_cbe_real_t qni, mni, wni;
  particle_cbe_real_t cni;
  // cni is just padding now. If I can managed to cut
  // one more element (possibly using the q/m ratio
  // modification) we can get rid of the padding all together. 
} particle_cbe_t;

#else

typedef struct psc_particle_cbe {
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

typedef struct psc_particles_cbe {
  particle_cbe_t *particles;
  int n_part;
} particles_cbe_t;

void particles_cbe_alloc(struct psc_particles *prts, particles_cbe_t *pp, int n_part);
void particles_cbe_realloc(struct psc_particles *prts, particles_cbe_t *pp, int new_n_part);
void particles_cbe_free(struct psc_particles *prts, particles_cbe_t *pp);

static inline particle_cbe_t *
particles_cbe_get_one(particles_cbe_t *pp, int n)
{
  return &pp->particles[n];
}

#endif

/// \file psc_particles_cbe.h Particle structure definitions and handeling
/// functions for the cell_be implementation. 
///
/// As Kai begins to implement his 'patch' system from libmrc, it becomes
/// less important that the cell have separate data structures. As this point, 
/// there are two things which prevent us from just being able to use the 
/// the default C particles base: we need the particles to be stored in 
/// 16 byte aligned memory; it makes life a great deal easier if each particle
/// is a multiple of 16 bytes. The memory alignment issue could easily be 
/// resolved by using posix_memalign in the C particles base. It might end up
/// being slower to allocate the memory (default is 8 in glibc), but it might
/// (possibly but unlikely) result in an overall speed increase. The issue 
/// of padding is significantly more difficult. Basically, we need to get rid
/// of one of the standard particle attributes. Obviously, the positions and 
/// momenta need to stay. That leaves us with somehow reducing qni,mni,and wni
/// to just two variables. If these two issues can be solved, we can do away with 
/// the PARTICLES_CBE base.

