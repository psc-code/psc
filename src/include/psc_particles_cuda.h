
#ifndef PSC_PARTICLE_CUDA_H
#define PSC_PARTICLE_CUDA_H

#include "psc_particles_private.h"
#include "cuda_wrap.h"
#include "psc_particles_single.h"

#define PTYPE PTYPE_CUDA
#include "psc_particles_common.h"
#undef PTYPE

struct cuda_bnd {
  psc_particle_cuda_buf_t buf;
  int n_recv;
  int n_send;
};

static inline int
particle_cuda_real_nint(particle_cuda_real_t x)
{
  return (int)(x + 10.5f) - 10;
}

static inline int
particle_cuda_real_fint(particle_cuda_real_t x)
{
  return (int)(x + 10.f) - 10;
}

#endif
