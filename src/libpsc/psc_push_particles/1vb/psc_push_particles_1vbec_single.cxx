
#include "psc_push_particles_private.h"

#include "push_particles_1vbec_single.hxx"

// ======================================================================
// psc_push_particles: subclass "1vbec_single"

using PushParticlesWrapper_t = PushParticlesWrapper<PushParticles1vbecSingle>;

struct psc_push_particles_ops_1vbec_single : psc_push_particles_ops {
  psc_push_particles_ops_1vbec_single() {
    name                  = "1vbec_single";
    size                  = PushParticlesWrapper_t::size;
    setup                 = PushParticlesWrapper_t::setup;
    destroy               = PushParticlesWrapper_t::destroy;
    particles_type        = "single";
  }
} psc_push_particles_1vbec_single_ops;

