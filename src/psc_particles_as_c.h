
#ifndef PSC_PARTICLES_AS_C_H
#define PSC_PARTICLES_AS_C_H

#include "psc_particles_c.h"

typedef particle_c_t particle_t;
typedef particles_c_t particles_t;

#define particles_get          particles_c_get
#define particles_put          particles_c_put
#define particles_get_one      particles_c_get_one

#endif

