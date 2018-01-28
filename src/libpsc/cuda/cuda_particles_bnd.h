
#ifndef CUDA_PARTICLES_BND_H
#define CUDA_PARTICLES_BND_H

#include "psc_particles_cuda.h"
#include "ddc_particles.hxx"

struct cuda_particles_bnd
{
  using mparticles_t = mparticles_cuda_t;
  using ddcp_t = ddc_particles<mparticles_t>;
  using ddcp_patch = typename ddcp_t::patch;

  void prep(ddcp_t* ddcp, cuda_mparticles* cmprts);
  void post(ddcp_t* ddcp, cuda_mparticles* cmprts);

  void spine_reduce(cuda_mparticles *cmprts);
  void spine_reduce_gold(cuda_mparticles *cmprts);
};

#endif

