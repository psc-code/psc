
#include "psc_push_particles_private.h"

#include "push_particles.hxx"
#include "push_config.hxx"
#include "push_dispatch.hxx"

using PushParticlesWrapper_t = PushParticlesWrapper<PushParticlesGenericC>;

// ======================================================================
// psc_push_particles: subclass "generic_c"

struct psc_push_particles_ops_c : psc_push_particles_ops {
  psc_push_particles_ops_c() {
    name                  = "generic_c";
    size                  = PushParticlesWrapper_t::size;
    setup                 = PushParticlesWrapper_t::setup;
    destroy               = PushParticlesWrapper_t::destroy;
  }
} psc_push_particles_generic_c_ops;
