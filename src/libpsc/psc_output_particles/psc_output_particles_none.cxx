
#include "psc_output_particles_private.h"
#include "output_particles.hxx"

struct psc_output_particles_none : OutputParticlesParams, OutputParticlesBase
{
  psc_output_particles_none(const Grid_t& grid, const OutputParticlesParams& params)
    : OutputParticlesParams(params)
  {}

  void run(MparticlesBase& mprts_base) override
  {}
};

// ======================================================================
// psc_output_particles: subclass "none"

struct psc_output_particles_ops_none : psc_output_particles_ops {
  using Wrapper_t = OutputParticlesWrapper<psc_output_particles_none>;
  psc_output_particles_ops_none() {
    name                  = "none";
    size                  = Wrapper_t::size;
    setup                 = Wrapper_t::setup;
    destroy               = Wrapper_t::destroy;
  }
} psc_output_particles_none_ops;
