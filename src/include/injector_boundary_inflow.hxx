#pragma once

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
