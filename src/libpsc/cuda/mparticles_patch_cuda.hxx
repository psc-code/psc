
#pragma once

#include "const_accessor_simple.hxx"

// FIXME?  we have two pretty similar versions of ConstAccessorCuda here,
// and probably only one should survive.  one version copies from
// device "on demand", whereas the other one copies a whole patch
// worth of data just in case.  Obviously, there can be use cases for
// either, but this needs some thinking as to what's actually needed
// for this code.  The on-demand version might be useful to serve as a
// template for a modifiable accesser, if we ever want to go there.

// ======================================================================
// ConstAccessorCuda

template <typename _Mparticles>
struct ConstAccessorCuda
{
  using Mparticles = _Mparticles;
  using _Particle = typename Mparticles::Particle;
  using Patch = ConstAccessorPatchSimple<ConstAccessorCuda>;

  ConstAccessorCuda(Mparticles& mprts)
    : mprts_{mprts},
      data_{const_cast<Mparticles&>(mprts).get_particles()},
      off_{mprts.get_offsets()}
  {}

  Patch operator[](int p) const { return {mprts_, p}; }
  Mparticles& mprts() const { return mprts_; }
  const _Particle* data(int p) const { return &data_[off_[p]]; }
  uint size(int p) const { return off_[p + 1] - off_[p]; }

private:
  Mparticles& mprts_;
  const std::vector<_Particle> data_;
  const std::vector<uint> off_;
};

// ======================================================================
// ConstPatchCuda

template <typename _Mparticles>
struct ConstPatchCuda
{
  using Mparticles = _Mparticles;

  ConstPatchCuda(const Mparticles& mprts, int p) : mprts_{mprts}, p_(p) {}

private:
  const Mparticles& mprts_;
  int p_;
};
