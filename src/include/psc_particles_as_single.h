
#ifndef PSC_PARTICLES_AS_SINGLE_H
#define PSC_PARTICLES_AS_SINGLE_H

#include "psc_particles_single.h"

typedef particle_single_real_t particle_real_t;
typedef particle_single_t particle_t;

#define particles_get_one           particles_single_get_one
#define particles_realloc           psc_particles_single_realloc
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

#define MPI_PARTICLES_REAL          MPI_PARTICLES_SINGLE_REAL
#define PARTICLE_TYPE               "single"

#define PSC_PARTICLES_AS_SINGLE 1

#include "particle_iter.h"

#endif

