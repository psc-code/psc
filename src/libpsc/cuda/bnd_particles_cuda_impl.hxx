
#pragma once

#include "bnd_particles_impl.hxx"
#include "psc_particles_cuda.h"

template<typename CudaMparticles, typename DIM>
struct cuda_bndp;

template<typename BS, typename DIM>
struct BndParticlesCuda : BndParticlesCommon<MparticlesCuda<BS>>
{
  using Base = BndParticlesCommon<MparticlesCuda<BS>>;

  BndParticlesCuda(struct mrc_domain *domain, const Grid_t& grid);
  ~BndParticlesCuda();

  void reset();
  void operator()(MparticlesCuda<BS>& mprts);
  void exchange_particles(MparticlesBase& mprts_base) override;

private:
  cuda_bndp<cuda_mparticles<BS>, DIM>* cbndp_;
};

