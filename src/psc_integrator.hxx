
#ifndef PSC_INTEGRATOR_HXX
#define PSC_INTEGRATOR_HXX

#include "psc.hxx"

namespace
{

class InjectParticlesNone
{
public:
  template <typename Mparticles>
  void operator()(const Grid_t& grid, Mparticles& mprts)
  {}
};

InjectParticlesNone injectParticlesNone;

} // namespace

// ======================================================================
// PscIntegrator

template <typename PscConfig, typename Diagnostics, typename InjectParticles>
struct PscIntegrator : Psc<PscConfig, Diagnostics>
{
  using Base = Psc<PscConfig, Diagnostics>;

  using MfieldsState = typename PscConfig::MfieldsState;
  using Mparticles = typename PscConfig::Mparticles;
  using Balance = typename PscConfig::Balance;
  using Collision = typename PscConfig::Collision;
  using Checks = typename PscConfig::Checks;
  using Marder = typename PscConfig::Marder;

  // ----------------------------------------------------------------------
  // ctor

  PscIntegrator(const PscParams& params, Grid_t& grid, MfieldsState& mflds,
                Mparticles& mprts, Balance& balance, Collision& collision,
                Checks& checks, Marder& marder, Diagnostics& diagnostics,
                InjectParticles& inject_particles)
    : Base{params,    grid,   mflds,  mprts,      balance,
           collision, checks, marder, diagnostics},
      inject_particles_{inject_particles}
  {}

  // ----------------------------------------------------------------------
  // inject_particles

  void inject_particles() override
  {
    return inject_particles_(Base::grid(), Base::mprts_);
  }

private:
  InjectParticles& inject_particles_;
};

template <typename PscConfig, typename MfieldsState, typename Mparticles,
          typename Balance, typename Collision, typename Checks,
          typename Marder, typename Diagnostics,
          typename InjectParticles = InjectParticlesNone>
PscIntegrator<PscConfig, Diagnostics, InjectParticles> makePscIntegrator(
  const PscParams& params, Grid_t& grid, MfieldsState& mflds, Mparticles& mprts,
  Balance& balance, Collision& collision, Checks& checks, Marder& marder,
  Diagnostics& diagnostics, InjectParticles& inject_particles = injectParticlesNone)
{
  return {params,    grid,   mflds,  mprts,       balance,
          collision, checks, marder, diagnostics, inject_particles};
}

#endif
