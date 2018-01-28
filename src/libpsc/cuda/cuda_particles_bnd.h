
#ifndef CUDA_PARTICLES_BND_H
#define CUDA_PARTICLES_BND_H

struct cuda_particles_bnd
{
  using mparticles_t = mparticles_cuda_t;
  using ddcp_t = ddc_particles<mparticles_t>;
  using ddcp_patch = typename ddcp_t::patch;

  void prep(ddcp_t* ddcp, cuda_mparticles* cmprts);
  void post(ddcp_t* ddcp, cuda_mparticles* cmprts);
};

#endif

