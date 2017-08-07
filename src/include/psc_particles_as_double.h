
#ifndef PSC_PARTICLES_AS_DOUBLE_H
#define PSC_PARTICLES_AS_DOUBLE_H

#include "psc_particles_double.h"

typedef particle_double_real_t particle_real_t;
typedef particle_double_t particle_t;

#define mparticles_get_one          psc_mparticles_double_get_one
#define particle_qni_div_mni        particle_double_qni_div_mni
#define particle_qni_wni            particle_double_qni_wni
#define particle_qni                particle_double_qni
#define particle_mni                particle_double_mni
#define particle_wni                particle_double_wni
#define particle_kind               particle_double_kind
#define particle_x                  particle_double_x
#define particle_px                 particle_double_px
#define particle_get_relative_pos   particle_double_get_relative_pos
#define particle_real_nint          particle_double_real_nint
#define particle_real_fint          particle_double_real_fint
#define particle_real_abs           particle_double_real_abs
#define particle_real_sqrt          particle_double_real_sqrt

#define MPI_PARTICLES_REAL          MPI_PARTICLES_DOUBLE_REAL
#define PARTICLE_TYPE               "double"

#define PSC_PARTICLES_AS_DOUBLE 1

#include "particle_iter.h"

#endif

