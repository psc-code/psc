
#pragma once

#include "bnd_particles_impl.hxx"
#include "psc_particles_cuda.h"

template<typename BS>
struct cuda_bndp;

template<typename BS>
struct BndParticlesCuda : BndParticlesCommon<MparticlesCuda>
{
  using Base = BndParticlesCommon<MparticlesCuda>;

  BndParticlesCuda(struct mrc_domain *domain, const Grid_t& grid);
  ~BndParticlesCuda();

  void reset();
  void operator()(MparticlesCuda& mprts);
  void exchange_particles(MparticlesBase& mprts_base) override;

private:
  cuda_bndp<BS>* cbndp_;
};

