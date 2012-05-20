
#ifndef PSC_PARTICLES_AS_FORTRAN_H
#define PSC_PARTICLES_AS_FORTRAN_H

#include "psc_particles_fortran.h"

typedef particle_fortran_t particle_t;
typedef mparticles_fortran_t mparticles_t;

#define psc_mparticles_get_from       psc_mparticles_fortran_get_from
#define psc_mparticles_put_to         psc_mparticles_fortran_put_to
#define particles_get_one      particles_fortran_get_one

#endif

