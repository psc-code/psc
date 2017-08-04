
#ifndef PSC_PARTICLES_AS_SINGLE_BY_BLOCK_H
#define PSC_PARTICLES_AS_SINGLE_BY_BLOCK_H

#include "psc_particles_single_by_block.h"

typedef particle_single_by_block_real_t particle_real_t;
typedef particle_single_by_block_t particle_t;

#define particles_get_one           particles_single_by_block_get_one
#define particles_realloc           particles_single_by_block_realloc
#define particle_qni_div_mni        particle_single_by_block_qni_div_mni
#define particle_qni_wni            particle_single_by_block_qni_wni
#define particle_qni                particle_single_by_block_qni
#define particle_mni                particle_single_by_block_mni
#define particle_wni                particle_single_by_block_wni
#define particle_kind               particle_single_by_block_kind
#define particle_x                  particle_single_by_block_x
#define particle_px                 particle_single_by_block_px
#define particle_get_relative_pos   particle_single_by_block_get_relative_pos
#define particle_real_nint          particle_single_by_block_real_nint
#define particle_real_fint          particle_single_by_block_real_fint
#define particle_real_abs           particle_single_by_block_real_abs
#define particle_real_sqrt          particle_single_by_block_real_sqrt

#define MPI_PARTICLES_REAL          MPI_PARTICLES_SINGLE_BY_BLOCK_REAL
#define PARTICLE_TYPE               "single_by_block"

#define PSC_PARTICLES_AS_SINGLE_BY_BLOCK 1

#include "particle_iter.h"

#endif

