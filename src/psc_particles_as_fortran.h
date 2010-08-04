
#ifndef PSC_PARTICLES_AS_FORTRAN_H
#define PSC_PARTICLES_AS_FORTRAN_H

#include "psc_particles_fortran.h"

typedef particle_fortran_t particle_t;
typedef psc_particles_fortran_t psc_particles_t;

#define psc_particles_get          psc_particles_fortran_get
#define psc_particles_put          psc_particles_fortran_put
#define psc_particles_get_one      psc_particles_fortran_get_one

#endif

