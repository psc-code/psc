
#ifndef PSC_PARTICLES_AS_CUDA2_H
#define PSC_PARTICLES_AS_CUDA2_H

#include "psc_particles_cuda2.h"

typedef particle_cuda2_real_t particle_real_t;
typedef particle_cuda2_t particle_t;

#define particles_get_one           particles_cuda2_get_one
#define particle_qni_div_mni        particle_cuda2_qni_div_mni
#define particle_qni_wni            particle_cuda2_qni_wni
#define particle_qni                particle_cuda2_qni
#define particle_mni                particle_cuda2_mni
#define particle_wni                particle_cuda2_wni
#define particle_kind               particle_cuda2_kind
#define particle_real_nint          particle_cuda2_real_nint

#define MPI_PARTICLES_REAL          MPI_PARTICLES_CUDA2_REAL
#define PARTICLE_TYPE               "cuda2"

#define PSC_PARTICLES_AS_CUDA2 1

#endif

