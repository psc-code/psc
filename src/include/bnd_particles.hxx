
#pragma once

#include "particles.hxx"

// ======================================================================
// BndParticlesBase

struct BndParticlesBase
{
  virtual void exchange_particles(MparticlesBase& mprts_base) = 0;
};

