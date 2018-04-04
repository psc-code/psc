
#pragma once

class CollisionCuda: CollisionBase
{
public:
  virtual void run(PscMparticlesBase mprts_base) { assert(0); }

  void operator()(MparticlesCuda& mprts)
  {
    MHERE;
  }
};

