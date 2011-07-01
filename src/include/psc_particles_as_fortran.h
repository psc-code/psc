
#ifndef PSC_PARTICLES_AS_FORTRAN_H
#define PSC_PARTICLES_AS_FORTRAN_H

#include "psc_particles_fortran.h"

typedef particle_fortran_t particle_t;
typedef particles_fortran_t particles_t;
typedef mparticles_fortran_t mparticles_t;

#define mparticles_get         mparticles_fortran_get
#define mparticles_put         mparticles_fortran_put
#define particles_get_one      particles_fortran_get_one

#endif

