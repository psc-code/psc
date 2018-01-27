
#include "psc_bnd_particles_private.h"
#include "psc_particles_single.h"

#include "psc_bnd_particles_common.cxx"

// ======================================================================
// psc_bnd_particles: subclass "single"

struct psc_bnd_particles_ops_single : psc_bnd_particles_ops {
  using sub = psc_bnd_particles_sub<mparticles_single_t>;
  psc_bnd_particles_ops_single() {
    name                    = "single";
    size                    = sizeof(sub);
    setup                   = sub::setup;
    unsetup                 = sub::unsetup;
    exchange_particles      = sub::exchange_particles;
    //open_calc_moments       = psc_bnd_particles_sub_open_calc_moments;
  }
} psc_bnd_particles_single_ops;
