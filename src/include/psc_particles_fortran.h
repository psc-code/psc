
#ifndef PSC_PARTICLE_FORTRAN_H
#define PSC_PARTICLE_FORTRAN_H

#include "psc_particles_private.h"

// this matches the Fortran particle data structure

typedef double particle_fortran_real_t;

#define MPI_PARTICLES_FORTRAN_REAL MPI_DOUBLE

typedef struct psc_particle_fortran {
  particle_fortran_real_t xi, yi, zi;
  particle_fortran_real_t pxi, pyi, pzi;
  particle_fortran_real_t qni;
  particle_fortran_real_t mni;
  particle_fortran_real_t cni;
  particle_fortran_real_t lni;
  particle_fortran_real_t wni;
} particle_fortran_t;

struct psc_particles_fortran {
  struct psc_particle_fortran *particles;
};

#define psc_particles_fortran(prts) mrc_to_subobj(prts, struct psc_particles_fortran)

static inline particle_fortran_t *
particles_fortran_get_one(struct psc_particles *prts, int n)
{
  assert(psc_particles_ops(prts) == &psc_particles_fortran_ops);
  return &psc_particles_fortran(prts)->particles[n];
}

static inline int
particle_fortran_real_fint(particle_fortran_real_t x)
{
  return (int)(x + 10.f) - 10;
}

#endif
