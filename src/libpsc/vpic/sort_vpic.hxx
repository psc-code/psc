
#pragma once

#include "psc_particles_vpic.h"

struct SortVpic
{
  using Mparticles = MparticlesVpic;

  void operator()(Mparticles& mprts)
  {
    Simulation_sort_mprts(mprts.sim, mprts.vmprts, ppsc->timestep);
  }
};

