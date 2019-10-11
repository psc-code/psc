
#pragma once

template<typename BS>
class SortCuda : SortBase
{
public:
  void operator()(MparticlesCuda<BS>& _mprts)
  {
    // nothing to be done, since MparticlesCuda are kept sorted anyway...
  }
};

