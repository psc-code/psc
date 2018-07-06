
#pragma once

template<typename BS>
class SortCuda : SortBase
{
public:
  virtual void run(MparticlesBase& mprts_base) { assert(0); }

  void operator()(MparticlesCuda<BS>& _mprts)
  {
    // nothing to be done, since MparticlesCuda are kept sorted anyway...
  }
};

