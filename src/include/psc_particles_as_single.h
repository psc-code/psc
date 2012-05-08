
#ifndef PSC_PARTICLES_AS_SINGLE_H
#define PSC_PARTICLES_AS_SINGLE_H

#include "psc_particles_single.h"

typedef particle_single_real_t particle_real_t;
typedef particle_single_t particle_t;
typedef particles_single_t particles_t;
typedef mparticles_single_t mparticles_t;

#define psc_mparticles_get_cf       psc_mparticles_get_single
#define psc_mparticles_put_cf       psc_mparticles_put_single
#define psc_mparticles_get_patch    psc_mparticles_get_patch_single
#define particles_get_one           particles_single_get_one
#define particles_realloc           particles_single_realloc
#define particle_qni_div_mni        particle_single_qni_div_mni
#define particle_qni_wni            particle_single_qni_wni
#define particle_qni                particle_single_qni
#define particle_wni                particle_single_wni
#define particle_get_relative_pos   particle_single_get_relative_pos
#define particle_real_nint          particle_single_real_nint
#define particle_real_fint          particle_single_real_fint

#define MPI_PARTICLES_REAL          MPI_PARTICLES_SINGLE_REAL

#endif

