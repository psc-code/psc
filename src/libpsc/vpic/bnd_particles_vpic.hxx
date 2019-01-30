
#pragma once

#include "bnd_particles.hxx"

template<typename Mparticles>
struct BndParticlesVpic : BndParticlesBase
{
  BndParticlesVpic(const Grid_t& grid) {}

  void operator()(Mparticles& mprts) {}
  void exchange_particles(MparticlesBase& mprts_base) override {}
};

