
#ifndef PSC_PARTICLES_AS_C_H
#define PSC_PARTICLES_AS_C_H

#include "psc_particles_c.h"

typedef particle_c_t particle_t;
typedef particles_c_t particles_t;
typedef mparticles_c_t mparticles_t;

#define mparticles_get         mparticles_c_get
#define mparticles_put         mparticles_c_put
#define particles_get_one      particles_c_get_one

#endif

