
#pragma once

#include "bnd_particles.hxx"
#include "psc_particles_vpic.h"

struct BndParticlesVpic : BndParticlesBase
{
  BndParticlesVpic(struct mrc_domain *domain, const Grid_t& grid) {}

  void operator()(MparticlesVpic& mprts) {}
  void exchange_particles(MparticlesBase& mprts_base) override {}
};

