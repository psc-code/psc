
#ifndef PSC_PARTICLES_AS_SINGLE_BY_BLOCK_H
#define PSC_PARTICLES_AS_SINGLE_BY_BLOCK_H

#include "psc_particles_single_by_block.h"

typedef particle_single_by_block_real_t particle_real_t;
typedef particle_single_by_block_t particle_t;

#define mparticles_get_one          psc_mparticles_single_by_block_get_one
#define mparticles_patch_get_b_mx   psc_mparticles_single_by_block_patch_get_b_mx
#define mparticles_patch_get_b_dxi  psc_mparticles_single_by_block_patch_get_b_dxi

#define particle_qni_div_mni        particle_single_by_block_qni_div_mni
#define particle_qni_wni            particle_single_by_block_qni_wni
#define particle_qni                particle_single_by_block_qni
#define particle_mni                particle_single_by_block_mni
#define particle_wni                particle_single_by_block_wni
#define particle_kind               particle_single_by_block_kind
#define particle_real_nint          particle_single_by_block_real_nint
#define particle_real_fint          particle_single_by_block_real_fint

#define particle_iter_t             psc_particle_single_by_block_iter_t
#define particle_iter_equal         psc_particle_single_by_block_iter_equal
#define particle_iter_next          psc_particle_single_by_block_iter_next
#define particle_iter_deref         psc_particle_single_by_block_iter_deref
#define particle_iter_at            psc_particle_single_by_block_iter_at
#define particle_range_t            psc_particle_single_by_block_range_t
#define particle_range_mprts        psc_particle_single_by_block_range_mprts
#define particle_range_size         psc_particle_single_by_block_range_size

#define MPI_PARTICLES_REAL          MPI_PARTICLES_SINGLE_BY_BLOCK_REAL
#define PARTICLE_TYPE               "single_by_block"

#define PSC_PARTICLES_AS_SINGLE_BY_BLOCK 1

#include "particle_iter.h"

#endif

