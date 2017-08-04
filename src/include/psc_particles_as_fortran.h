
#ifndef PSC_PARTICLES_AS_FORTRAN_H
#define PSC_PARTICLES_AS_FORTRAN_H

#include "psc_particles_fortran.h"

typedef particle_fortran_real_t particle_real_t;
typedef particle_fortran_t particle_t;

#define particles_get_one             particles_fortran_get_one
#define particles_realloc             particles_fortran_realloc
#define particle_real_fint            particle_fortran_real_fint

#define MPI_PARTICLES_REAL            MPI_PARTICLES_FORTRAN_REAL
#define PARTICLE_TYPE                 "fortran"

#include "particle_iter.h"

#endif

