
#include "output_particles_none_impl.hxx"

// ======================================================================
// psc_output_particles: subclass "none"

struct psc_output_particles_ops_none : psc_output_particles_ops {
  using Wrapper_t = OutputParticlesWrapper<OutputParticlesNone>;
  psc_output_particles_ops_none() {
    name                  = "none";
    size                  = Wrapper_t::size;
    setup                 = Wrapper_t::setup;
    destroy               = Wrapper_t::destroy;
  }
} psc_output_particles_none_ops;
