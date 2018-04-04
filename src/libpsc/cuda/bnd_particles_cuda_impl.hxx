
#pragma once

#include "bnd_particles_impl.hxx"
#include "psc_particles_cuda.h"

struct cuda_bndp;

struct BndParticlesCuda : BndParticles_<MparticlesCuda>
{
  using Base = BndParticles_<MparticlesCuda>;

  BndParticlesCuda(struct mrc_domain *domain, const Grid_t& grid);
  ~BndParticlesCuda();

  void reset();
  void exchange_particles(PscMparticlesBase mprts_base) override;

private:
  cuda_bndp* cbndp_;
};

