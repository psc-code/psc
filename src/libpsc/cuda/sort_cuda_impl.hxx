
#pragma once

class SortCuda : SortBase
{
public:
  virtual void run(PscMparticlesBase mprts_base) { assert(0); }

  void operator()(MparticlesCuda<BS144>& _mprts)
  {
    // nothing to be done, since MparticlesCuda are kept sorted anyway...
  }
};

