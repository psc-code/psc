
#ifndef PSC_PARTICLES_AS_SINGLE_BY_BLOCK_H
#define PSC_PARTICLES_AS_SINGLE_BY_BLOCK_H

#include "psc_particles_single_by_block.h"

typedef particle_single_by_block_real_t particle_real_t;
typedef particle_single_by_block_t particle_t;
using mparticles_t = mparticles_single_by_block_t;

#define particle_qni_div_mni        particle_single_by_block_qni_div_mni
#define particle_qni_wni            particle_single_by_block_qni_wni
#define particle_qni                particle_single_by_block_qni
#define particle_mni                particle_single_by_block_mni
#define particle_wni                particle_single_by_block_wni

#define particle_iter_t             psc_particle_single_by_block_iter_t
#define particle_range_t            psc_particle_single_by_block_range_t

#define PARTICLE_TYPE               "single_by_block"

#define PSC_PARTICLES_AS_SINGLE_BY_BLOCK 1

#include "particle_iter.h"

#endif

