
#pragma once

#include "bnd_particles.hxx"
#include "vpic_iface.h"

struct BndParticlesVpic : BndParticlesBase
{
  BndParticlesVpic(struct mrc_domain *domain, const Grid_t& grid) {}

  void exchange_particles(MparticlesBase& mprts_base) override {}
};

