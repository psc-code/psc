
#ifndef PSC_PARTICLES_AS_SINGLE_H
#define PSC_PARTICLES_AS_SINGLE_H

#include "psc_particles_single.h"

typedef particle_single_real_t particle_real_t;
typedef particle_single_t particle_t;

#define mparticles_get_one          psc_mparticles_single_get_one
#define mparticles_get_n_prts       psc_mparticles_single_get_n_prts
#define mparticles_push_back        psc_mparticles_single_push_back

#define particle_qni_div_mni        particle_single_qni_div_mni
#define particle_qni_wni            particle_single_qni_wni
#define particle_qni                particle_single_qni
#define particle_mni                particle_single_mni
#define particle_wni                particle_single_wni
#define particle_kind               particle_single_kind
#define particle_x                  particle_single_x
#define particle_px                 particle_single_px
#define particle_get_relative_pos   particle_single_get_relative_pos
#define particle_real_nint          particle_single_real_nint
#define particle_real_fint          particle_single_real_fint
#define particle_real_abs           particle_single_real_abs
#define particle_real_sqrt          particle_single_real_sqrt

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

