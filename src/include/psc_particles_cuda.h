
#ifndef PSC_PARTICLE_CUDA_H
#define PSC_PARTICLE_CUDA_H

#include "psc_particles_private.h"
#include "psc_particles_single.h"

#include "psc_particle_buf_cuda.h"

#define PTYPE PTYPE_CUDA
#include "psc_particles_common.h"
#undef PTYPE

static inline int
particle_cuda_real_nint(particle_cuda_real_t x)
{
  return (int)(x + 10.5f) - 10;
}

#endif
