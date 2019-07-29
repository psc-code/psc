
#ifndef PSC_INTEGRATOR_HXX
#define PSC_INTEGRATOR_HXX

namespace
{

template <typename Mparticles>
void injectNone(const Grid_t& grid, Mparticles& mprts)
{}

} // namespace

// ======================================================================
// PscIntegrator

template <typename PscConfig, typename InjectFunc>
struct PscIntegrator : Psc<PscConfig>
{
  using Base = Psc<PscConfig>;

  using MfieldsState = typename PscConfig::MfieldsState;
  using Mparticles = typename PscConfig::Mparticles;
  using Balance = typename PscConfig::Balance;
  using Collision = typename PscConfig::Collision;
  using Checks = typename PscConfig::Checks;
  using Marder = typename PscConfig::Marder;
  using OutputParticles = typename PscConfig::OutputParticles;
  using Diagnostics = typename Base::Diagnostics;

  // ----------------------------------------------------------------------
  // ctor

  PscIntegrator(const PscParams& params, Grid_t& grid, MfieldsState& mflds,
                Mparticles& mprts, Balance& balance, Collision& collision,
                Checks& checks, Marder& marder, Diagnostics& diagnostics,
                InjectFunc& inject_particles)
    : Base{params, grid,   mflds, mprts, balance,    collision,
           checks, marder, diagnostics},
      inject_particles_{inject_particles}
  {}

  // ----------------------------------------------------------------------
  // inject_particles

  void inject_particles() override
  {
    return inject_particles_(Base::grid(), Base::mprts_);
  }

private:
  InjectFunc& inject_particles_;
};

template <typename PscConfig, typename MfieldsState, typename Mparticles,
          typename Balance, typename Collision, typename Checks,
          typename Marder, typename Diagnostics, typename InjectFunc>
PscIntegrator<PscConfig, InjectFunc> makePscIntegrator(
  const PscParams& params, Grid_t& grid, MfieldsState& mflds, Mparticles& mprts,
  Balance& balance, Collision& collision, Checks& checks, Marder& marder,
  Diagnostics& diagnostics, InjectFunc& inject_particles)
{
  return {params,    grid,   mflds,  mprts,       balance,
          collision, checks, marder, diagnostics, inject_particles};
}

template <typename PscConfig, typename MfieldsState, typename Mparticles,
          typename Balance, typename Collision, typename Checks,
          typename Marder, typename Diagnostics>
PscIntegrator<PscConfig, decltype(injectNone<Mparticles>)> makePscIntegrator(
  const PscParams& params, Grid_t& grid, MfieldsState& mflds, Mparticles& mprts,
  Balance& balance, Collision& collision, Checks& checks, Marder& marder,
  Diagnostics& diagnostics)
{
  return {params,    grid,   mflds,  mprts,       balance,
          collision, checks, marder, diagnostics, injectNone<Mparticles>};
}

#endif
