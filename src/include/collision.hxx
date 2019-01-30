
#pragma once

#include "particles.hxx"

// ======================================================================
// CollisionBase

class CollisionBase
{
public:
  virtual void operator()(MparticlesBase& mprts_base) = 0;
};

