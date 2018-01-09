
#ifndef PSC_PARTICLES_CUDA_H
#define PSC_PARTICLES_CUDA_H

#include "psc_particles_private.h"
#include "psc_particles_single.h"

#include "psc_particle_buf_cuda.h"

#define PTYPE PTYPE_CUDA
#include "psc_particles_common.h"
#undef PTYPE

#endif
