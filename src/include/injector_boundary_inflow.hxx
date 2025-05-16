#pragma once

#include <math.h>

#include "grid.hxx"
#include "rng.hxx"
#include "particle.h"

class ParticleGeneratorMaxwellian
{
public:
  using Real = psc::particle::Inject::Real;
  using Real3 = psc::particle::Inject::Real3;

  // FIXME would be nice to just pass 1 thing for kind-related info
  ParticleGeneratorMaxwellian(int kind_idx, Grid_t::Kind kind, Real w,
                              Real3 mean_u, Real3 temperature)
    : kind_idx{kind_idx}, w{w}
  {
    for (int d = 0; d < 3; d++) {
      Real stdev_u = sqrt(temperature[d] / kind.m);
      vdfs[d] = VelocityDistributionFunction{mean_u[d], stdev_u};
    }
  }

  psc::particle::Inject get()
  {
    // TODO: set x later, since we don't have any grid info here

    Real3 x{0.0, 0.0, 0.0};
    Real3 u{vdfs[0].get(), vdfs[1].get(), vdfs[2].get()};
    psc::particle::Tag tag = 0;

    return {x, u, w, kind_idx, tag};
  }

private:
  using VelocityDistributionFunction = rng::Normal<Real>;
  Vec3<VelocityDistributionFunction> vdfs;
  Real w;
  int kind_idx;
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
