
#ifndef PSC_PARTICLES_AS_C_H
#define PSC_PARTICLES_AS_C_H

#include "psc_particles_c.h"

typedef particle_c_real_t particle_real_t;
typedef particle_c_t particle_t;

#define mparticles_get_one          psc_mparticles_c_get_one
#define mparticles_get_n_prts       psc_mparticles_c_get_n_prts
#define mparticles_patch_reserve    psc_mparticles_c_patch_reserve
#define mparticles_patch_push_back  psc_mparticles_c_patch_push_back

#define particle_qni_div_mni        particle_c_qni_div_mni
#define particle_qni_wni            particle_c_qni_wni
#define particle_qni                particle_c_qni
#define particle_mni                particle_c_mni
#define particle_wni                particle_c_wni
#define particle_kind               particle_c_kind
#define particle_x                  particle_c_x
#define particle_px                 particle_c_px
#define particle_get_relative_pos   particle_c_get_relative_pos
#define particle_real_nint          particle_c_real_nint
#define particle_real_fint          particle_c_real_fint
#define particle_real_sqrt          particle_c_real_sqrt
#define particle_real_abs           particle_c_real_abs

#define particle_iter_t             psc_particle_c_iter_t
#define particle_iter_equal         psc_particle_c_iter_equal
#define particle_iter_next          psc_particle_c_iter_next 
#define particle_iter_deref         psc_particle_c_iter_deref 
#define particle_iter_at            psc_particle_c_iter_at
#define particle_range_t            psc_particle_c_range_t
#define particle_range_mprts        psc_particle_c_range_mprts
#define particle_range_size         psc_particle_c_range_size

#define MPI_PARTICLES_REAL          MPI_PARTICLES_C_REAL
#define PARTICLE_TYPE               "c"

#define PSC_PARTICLES_AS_C 1

#include "particle_iter.h"

#endif

