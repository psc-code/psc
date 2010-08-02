
#ifndef PSC_PARTICLE_FORTRAN_H
#define PSC_PARTICLE_FORTRAN_H

#include "psc.h"

// this matches the Fortran particle data structure

typedef f_real particle_fortran_real_t;


typedef struct {
  particle_fortran_real_t xi, yi, zi;
  particle_fortran_real_t pxi, pyi, pzi;
  particle_fortran_real_t qni;
  particle_fortran_real_t mni;
  particle_fortran_real_t cni;
  particle_fortran_real_t lni;
  particle_fortran_real_t wni;
} particle_fortran_t;

#endif
