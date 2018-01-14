
#ifndef PSC_PARTICLES_AS_SINGLE_H
#define PSC_PARTICLES_AS_SINGLE_H

#include "psc_particles_single.h"

typedef particle_single_real_t particle_real_t;
using mparticles_t = mparticles_single_t;
using particle_t = mparticles_t::particle_t;

#define particle_qni_div_mni        particle_single_qni_div_mni
#define particle_qni_wni            particle_single_qni_wni
#define particle_qni                particle_single_qni
#define particle_mni                particle_single_mni
#define particle_wni                particle_single_wni

#define PARTICLE_TYPE               "single"

#define PSC_PARTICLES_AS_SINGLE 1

#include "particle_iter.h"

#endif

