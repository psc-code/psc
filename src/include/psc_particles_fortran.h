
#ifndef PSC_PARTICLE_FORTRAN_H
#define PSC_PARTICLE_FORTRAN_H

#include "psc.h"

// this matches the Fortran particle data structure

typedef double particle_fortran_real_t;

#define MPI_PARTICLES_FORTRAN_REAL MPI_DOUBLE

typedef struct {
  particle_fortran_real_t xi, yi, zi;
  particle_fortran_real_t pxi, pyi, pzi;
  particle_fortran_real_t qni;
  particle_fortran_real_t mni;
  particle_fortran_real_t cni;
  particle_fortran_real_t lni;
  particle_fortran_real_t wni;
} particle_fortran_t;

typedef struct {
  particle_fortran_t *particles;
  int n_part;
  int n_alloced;
} particles_fortran_t;

void particles_fortran_alloc(particles_fortran_t *pp, int n_part);
void particles_fortran_realloc(particles_fortran_t *pp, int new_n_part);
void particles_fortran_free(particles_fortran_t *pp);

static inline particle_fortran_t *
particles_fortran_get_one(particles_fortran_t *pp, int n)
{
  return &pp->particles[n];
}

#endif
