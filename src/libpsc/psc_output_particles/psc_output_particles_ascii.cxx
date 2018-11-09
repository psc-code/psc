
#include "output_particles_ascii_impl.hxx"

// ======================================================================
// psc_output_particles: subclass "ascii"

struct psc_output_particles_ops_ascii : psc_output_particles_ops {
  using Wrapper_t = OutputParticlesWrapper<psc_output_particles_ascii>;
  psc_output_particles_ops_ascii() {
    name                  = "ascii";
    size                  = Wrapper_t::size;
    setup                 = Wrapper_t::setup;
    destroy               = Wrapper_t::destroy;
  }
} psc_output_particles_ascii_ops;
