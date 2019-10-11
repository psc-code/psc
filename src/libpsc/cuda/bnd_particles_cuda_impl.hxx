
#pragma once

#include "bnd_particles_impl.hxx"

template<typename CudaMparticles, typename DIM>
struct cuda_bndp;

template<typename Mparticles, typename DIM>
struct BndParticlesCuda : BndParticlesCommon<Mparticles>
{
  using Base = BndParticlesCommon<Mparticles>;

  BndParticlesCuda(const Grid_t& grid);
  ~BndParticlesCuda();

  void reset(const Grid_t& grid);
  void operator()(Mparticles& mprts);

private:
  cuda_bndp<typename Mparticles::CudaMparticles, DIM>* cbndp_;
};

