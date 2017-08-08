
#ifndef PSC_PARTICLES_AS_FORTRAN_H
#define PSC_PARTICLES_AS_FORTRAN_H

#include "psc_particles_fortran.h"

typedef particle_fortran_real_t particle_real_t;
typedef particle_fortran_t particle_t;

#define mparticles_get_one            psc_mparticles_fortran_get_one
#define particle_real_fint            particle_fortran_real_fint

#define particle_iter_t               psc_particle_fortran_iter_t

#define MPI_PARTICLES_REAL            MPI_PARTICLES_FORTRAN_REAL
#define PARTICLE_TYPE                 "fortran"

#define PSC_PARTICLES_AS_FORTRAN 1

#include "particle_iter.h"

#endif

