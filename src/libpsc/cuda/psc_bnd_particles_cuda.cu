
#include "psc_bnd_particles_private.h"
#include "psc_particles_cuda.h"

#include "cuda_bndp_iface.h"

#include <mrc_profile.h>

using mparticles_t = mparticles_cuda_t;

// ======================================================================
// psc_bnd_particles: subclass "cuda"

struct psc_bnd_particles_ops_cuda : psc_bnd_particles_ops {
  psc_bnd_particles_ops_cuda() {
    name                    = "cuda";
    size                    = sizeof(psc_bnd_particles_cuda);
    destroy                 = psc_bnd_particles_cuda::destroy;
    setup                   = psc_bnd_particles_cuda::setup;
    reset                   = psc_bnd_particles_cuda::reset;
    exchange_particles      = psc_bnd_particles_cuda::exchange_particles;
  }
} psc_bnd_particles_cuda_ops;

