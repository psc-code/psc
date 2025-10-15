#pragma once

#include "grid.hxx"
#include "particle.h"

class ParticleGeneratorMaxwellian
{
  using Real = psc::particle::Inject::Real;
  using Real3 = psc::particle::Inject::Real3;

public:
  // FIXME would be nice to just pass 1 thing for kind-related info
  ParticleGeneratorMaxwellian(int kind_idx, Grid_t::Kind kind, Real3 mean_u,
                              Real3 temperature)
  {}

  psc::particle::Inject get()
  {

    // TODO:
    // 1. sample x, u
    // 2. get w, kind, tag

    Real3 x{0.0, 0.0, 0.0};
    Real3 u{0.0, 0.0, 0.0};
    Real w{0.0};
    int kind = 0;
    psc::particle::Tag tag = 0;

    return {x, u, w, kind, tag};
  }
};

class InjectorBoundaryInflow
{
public:
  template <typename Mparticles, typename MfieldsState>
  void operator()(Mparticles& mprts, MfieldsState& mflds)
  {
    // TODO:
    // 1. inject particles
    // 2. advance injected particles
    // 3. cull injected particles that don't enter the domain
    // 4. update current
  }
};
