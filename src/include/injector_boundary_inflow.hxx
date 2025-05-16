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

  psc::particle::Inject get(Real3 min_pos, Real3 pos_range)
  {
    Real3 x;
    for (int d = 0; d < 3; d++) {
      x[d] = min_pos[d] + uniform_dist.get() * pos_range[d];
    }

    Real3 u{vdfs[0].get(), vdfs[1].get(), vdfs[2].get()};
    psc::particle::Tag tag = 0;

    return {x, u, w, kind_idx, tag};
  }

private:
  using VelocityDistributionFunction = rng::Normal<Real>;
  Vec3<VelocityDistributionFunction> vdfs;
  Real w;
  int kind_idx;
  rng::Uniform<Real> uniform_dist{0.0, 1.0};
};

template <typename PARTICLE_GENERATOR, typename PUSH_PARTICLES>
class InjectorBoundaryInflow
{
public:
  using ParticleGenerator = PARTICLE_GENERATOR;
  using PushParticles = PUSH_PARTICLES;

  using Mparticles = typename PushParticles::Mparticles;
  using MfieldsState = typename PushParticles::MfieldsState;
  using AdvanceParticle_t = typename PushParticles::AdvanceParticle_t;
  using InterpolateEM_t = typename PushParticles::InterpolateEM_t;
  using Current = typename PushParticles::Current;
  using Dim = typename PushParticles::Dim;
  using Real = typename PushParticles::real_t;
  using Real3 = typename PushParticles::Real3;
  using checks_order = typename PushParticles::checks_order;

  InjectorBoundaryInflow(ParticleGenerator particle_generator)
    : partice_generator{particle_generator}
  {}

  void operator()(Mparticles& mprts, MfieldsState& mflds)
  {
    // TODO:
    // 1. inject particles
    // 2. advance injected particles
    // 3. cull injected particles that don't enter the domain
    // 4. update current
  }

private:
  ParticleGenerator partice_generator;
};
