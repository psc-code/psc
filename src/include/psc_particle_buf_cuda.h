
#ifndef PSC_PARTICLE_BUF_CUDA_H
#define PSC_PARTICLE_BUF_CUDA_H

#include "particles.hxx"

#include <vector>

using particle_cuda_real_t = float;

struct particle_cuda_t : psc_particle<particle_cuda_real_t> {};

using psc_particle_cuda_buf_t = std::vector<particle_cuda_t>;

#endif
