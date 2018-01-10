
#ifndef PSC_PARTICLE_BUF_CUDA_H
#define PSC_PARTICLE_BUF_CUDA_H

#include "particles.hxx"

using particle_cuda_real_t = float;

struct particle_cuda_t : psc_particle<particle_cuda_real_t> {};

using psc_particle_cuda_buf_t = psc_particle_buf<particle_cuda_t>;

#define psc_mparticles_cuda(mprts) mrc_to_subobj(mprts, struct psc_mparticles_cuda)

#endif
