
#ifndef PSC_PARTICLES_AS_SINGLE_H
#define PSC_PARTICLES_AS_SINGLE_H

#include "psc_particles_single.h"

typedef particle_single_real_t particle_real_t;
typedef particle_single_t particle_t;

#define mparticles_patch_get_buf    psc_mparticles_single_patch_get_buf
#define mparticles_get_one          psc_mparticles_single_get_one
#define mparticles_get_n_prts       psc_mparticles_single_get_n_prts
#define mparticles_patch_reserve    psc_mparticles_single_patch_reserve
#define mparticles_patch_push_back  psc_mparticles_single_patch_push_back
#define mparticles_patch_resize     psc_mparticles_single_patch_resize
#define mparticles_patch_capacity   psc_mparticles_single_patch_capacity
#define mparticles_patch_get_b_mx   psc_mparticles_single_patch_get_b_mx
#define mparticles_patch_get_b_dxi  psc_mparticles_single_patch_get_b_dxi

#define particle_buf_t              psc_particle_single_buf_t
#define particle_buf_t              psc_particle_single_buf_t
#define particle_buf_ctor           psc_particle_single_buf_ctor
#define particle_buf_dtor           psc_particle_single_buf_dtor
#define particle_buf_size           psc_particle_single_buf_size
#define particle_buf_resize         psc_particle_single_buf_resize
#define particle_buf_reserve        psc_particle_single_buf_reserve
#define particle_buf_capacity       psc_particle_single_buf_capacity
#define particle_buf_push_back      psc_particle_single_buf_push_back
#define particle_buf_at_ptr         psc_particle_single_buf_at_ptr

#define particle_qni_div_mni        particle_single_qni_div_mni
#define particle_qni_wni            particle_single_qni_wni
#define particle_qni                particle_single_qni
#define particle_mni                particle_single_mni
#define particle_wni                particle_single_wni
#define particle_kind               particle_single_kind

#define particle_iter_t             psc_particle_single_iter_t
#define particle_iter_equal         psc_particle_single_iter_equal
#define particle_iter_next          psc_particle_single_iter_next
#define particle_iter_deref         psc_particle_single_iter_deref
#define particle_iter_at            psc_particle_single_iter_at
#define particle_range_t            psc_particle_single_range_t
#define particle_range_mprts        psc_particle_single_range_mprts
#define particle_range_size         psc_particle_single_range_size

#define MPI_PARTICLES_REAL          MPI_PARTICLES_SINGLE_REAL
#define PARTICLE_TYPE               "single"

#define PSC_PARTICLES_AS_SINGLE 1

#include "particle_iter.h"

#endif

