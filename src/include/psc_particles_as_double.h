
#ifndef PSC_PARTICLES_AS_DOUBLE_H
#define PSC_PARTICLES_AS_DOUBLE_H

#include "psc_particles_double.h"

typedef particle_double_real_t particle_real_t;
typedef particle_double_t particle_t;
using mparticles_t = mparticles_double_t;

#define particle_buf_t              psc_particle_double_buf_t

#define particle_qni_div_mni        particle_double_qni_div_mni
#define particle_qni_wni            particle_double_qni_wni
#define particle_qni                particle_double_qni
#define particle_mni                particle_double_mni
#define particle_wni                particle_double_wni

#define particle_iter_t             psc_particle_double_iter_t
#define particle_range_t            psc_particle_double_range_t

#define PARTICLE_TYPE               "double"

#define PSC_PARTICLES_AS_DOUBLE 1

#include "particle_iter.h"

#endif

