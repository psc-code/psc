
#include "psc_bnd_particles_private.h"
#include "psc_particles_as_single.h"

#include "psc_bnd_particles_common2.cxx"

// ======================================================================
// psc_bnd_particles: subclass "single2"

struct psc_bnd_particles_ops_single2 : psc_bnd_particles_ops {
  using sub_t = psc_bnd_particles_sub<mparticles_single_t,
				      bnd_particles_policy_ordered<mparticles_single_t>>;
  psc_bnd_particles_ops_single2() {
    name                    = "single2";
    size                    = sizeof(sub_t);
    setup                   = sub_t::setup;
    unsetup                 = sub_t::unsetup;
    exchange_particles      = sub_t::exchange_particles;
  }
} psc_bnd_particles_single2_ops;
