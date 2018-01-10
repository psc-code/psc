
#ifndef PSC_PARTICLE_BUF_CUDA_H
#define PSC_PARTICLE_BUF_CUDA_H

#include "particles.hxx"

#define PTYPE PTYPE_CUDA
#include "psc_particle_common.h"
#undef PTYPE

using psc_particle_cuda_buf_t = psc_particle_buf<particle_cuda_t>;

#endif
