
#ifndef PSC_PARTICLES_AS_C_H
#define PSC_PARTICLES_AS_C_H

#include "psc_particles_c.h"

typedef particle_c_real_t particle_real_t;
typedef particle_c_t particle_t;
typedef particles_c_t particles_t;
typedef mparticles_c_t mparticles_t;

#define psc_mparticles_create       psc_mparticles_c_create
#define psc_mparticles_set_domain_nr_particles psc_mparticles_c_set_domain_nr_particles
#define psc_mparticles_setup        psc_mparticles_c_setup
#define psc_mparticles_base_get_cf  psc_mparticles_base_get_c
#define psc_mparticles_base_put_cf  psc_mparticles_base_put_c
#define particles_get_one           particles_c_get_one
#define particles_realloc           particles_c_realloc
#define particle_real_nint          particle_c_real_nint
#define particle_real_fint          particle_c_real_fint

#define MPI_PARTICLES_REAL          MPI_PARTICLES_C_REAL

#endif

