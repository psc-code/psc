
#include "psc_bnd_particles_private.h"
#include "psc_particles_double.h"

#include "psc_bnd_particles_common.cxx"

// ======================================================================
// psc_bnd_particles: subclass "double"

struct psc_bnd_particles_ops_double : psc_bnd_particles_ops {
  using sub = psc_bnd_particles_sub<mparticles_double_t>;
  psc_bnd_particles_ops_double() {
    name                    = "double";
    size                    = sizeof(sub);
    setup                   = sub::setup;
    unsetup                 = sub::unsetup;
    exchange_particles      = sub::exchange_particles;
    //open_calc_moments       = psc_bnd_particles_sub_open_calc_moments;
  }
} psc_bnd_particles_double_ops;
