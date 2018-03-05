
#include "psc_bnd_particles_private.h"
#include "psc_particles_cuda.h"

#include "cuda_bndp_iface.h"

#include <mrc_profile.h>

using mparticles_t = PscMparticlesCuda;

// ======================================================================
// psc_bnd_particles: subclass "cuda"

struct psc_bnd_particles_ops_cuda : psc_bnd_particles_ops {
  using PscBndParticles_t = PscBndParticlesWrapper<psc_bnd_particles_cuda>;
  psc_bnd_particles_ops_cuda() {
    name                    = "cuda";
    size                    = sizeof(psc_bnd_particles_cuda);
    destroy                 = PscBndParticles_t::destroy;
    setup                   = PscBndParticles_t::setup;
    reset                   = PscBndParticles_t::reset;
    exchange_particles      = psc_bnd_particles_cuda::exchange_particles;
  }
} psc_bnd_particles_cuda_ops;

