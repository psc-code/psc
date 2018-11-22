
#pragma once

// ======================================================================
// PatchCuda

template<typename _Mparticles>
struct PatchCuda
{
  using Mparticles = _Mparticles;

  PatchCuda(Mparticles& mprts, int p)
    : mprts_(mprts), p_(p)
  {}
  
  const Grid_t& grid() const { return mprts_.grid(); }

protected:
  Mparticles& mprts_;
  int p_;
};

