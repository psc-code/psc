
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

struct psc_mparticles_fortran_patch {
  particle_fortran_t *prt_array;
};

struct psc_mparticles_fortran {
  struct psc_mparticles_fortran_patch *patch;
};

#define psc_mparticles_fortran(prts) mrc_to_subobj(prts, struct psc_mparticles_fortran)

static inline particle_fortran_t *
particles_fortran_get_one(struct psc_particles *prts, int n)
{
  struct psc_mparticles *mprts = prts->mprts;
  int p = prts->p;
  
  assert(psc_mparticles_ops(mprts) == &psc_mparticles_fortran_ops);
  return &psc_mparticles_fortran(mprts)->patch[p].prt_array[n];
}

static inline int
particle_fortran_real_fint(particle_fortran_real_t x)
{
  return (int)(x + 10.f) - 10;
}

#endif
