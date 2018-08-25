
#pragma once

#include "bnd_particles.hxx"

struct BndParticlesVpic : BndParticlesBase
{
  BndParticlesVpic(const Grid_t& grid) {}

  void operator()(MparticlesVpic& mprts) {}
  void exchange_particles(MparticlesBase& mprts_base) override {}
};

