
#ifndef PSC_INTEGRATOR_HXX
#define PSC_INTEGRATOR_HXX

// ======================================================================
// PscIntegrator

template <typename InjectFunc, typename PscConfig>
struct PscIntegrator : Psc<PscConfig>
{
  using Base = Psc<PscConfig>;

  using MfieldsState = typename PscConfig::MfieldsState;
  using Mparticles = typename PscConfig::Mparticles;
  using Balance = typename PscConfig::Balance_t;
  using Collision = typename PscConfig::Collision_t;
  using Checks = typename PscConfig::Checks_t;
  using Marder = typename PscConfig::Marder_t;
  using OutputParticles = typename PscConfig::OutputParticles;

  // ----------------------------------------------------------------------
  // ctor

  PscIntegrator(const PscParams& params, Grid_t& grid, MfieldsState& mflds,
                Mparticles& mprts, Balance& balance, Collision& collision,
                Checks& checks, Marder& marder, OutputFieldsC& outf,
                OutputParticles& outp, InjectFunc& inject_particles)
    : inject_particles_{inject_particles}
  {
    auto comm = grid.comm();

    Base::p_ = params;

    Base::define_grid(grid);
    Base::define_field_array(mflds);
    Base::define_particles(mprts);

    Base::balance_.reset(&balance);
    Base::collision_.reset(&collision);
    Base::checks_.reset(&checks);
    Base::marder_.reset(&marder);

    Base::outf_.reset(&outf);
    Base::outp_.reset(&outp);

    Base::init();

    Base::initialize();
  }

  // ----------------------------------------------------------------------
  // inject_particles

  void inject_particles() override
  {
    return inject_particles_(Base::grid(), *Base::mprts_);
  }

private:
  InjectFunc& inject_particles_;
};

#endif
