
#pragma once

#include "psc_particles_vpic.h"

struct SortVpic
{
  using Mparticles = MparticlesVpic;

  void operator()(Mparticles& mprts)
  {
    mprts.sim_->sort_mprts(mprts.vmprts_, ppsc->timestep);
  }
};

