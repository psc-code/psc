
#include "psc_bnd_particles_private.h"
#include "psc_particles_as_single.h"

#include "psc_bnd_particles_common2.cxx"

// ======================================================================
// psc_bnd_particles: subclass "single2"

struct psc_bnd_particles_ops_single2 : psc_bnd_particles_ops {
  psc_bnd_particles_ops_single2() {
    name                    = "single2";
    setup                   = psc_bnd_particles_sub_setup;
    unsetup                 = psc_bnd_particles_sub_unsetup;
    exchange_mprts_prep     = mparticles_ddcp<mparticles_t>::exchange_mprts_prep;
    exchange_mprts_post     = mparticles_ddcp<mparticles_t>::exchange_mprts_post;
    exchange_particles      = psc_bnd_particles_sub_exchange_particles;
  }
} psc_bnd_particles_single2_ops;
