
#include "psc_bnd_particles_private.h"
#include "psc_particles_as_double.h"

#include "psc_bnd_particles_common.cxx"

// ======================================================================
// psc_bnd_particles: subclass "double"

struct psc_bnd_particles_ops_double : psc_bnd_particles_ops {
  psc_bnd_particles_ops_double() {
    name                    = "double";
    setup                   = psc_bnd_particles_sub_setup;
    unsetup                 = psc_bnd_particles_sub_unsetup;
    exchange_particles      = psc_bnd_particles_sub_exchange_particles;
    open_calc_moments       = psc_bnd_particles_sub_open_calc_moments;
  }
} psc_bnd_particles_double_ops;
