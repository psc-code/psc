
#include "psc_bnd_particles_private.h"
#include "psc_particles_as_single.h"

#include "psc_bnd_particles_common.cxx"

// ======================================================================
// psc_bnd_particles: subclass "single"

struct psc_bnd_particles_ops_single : psc_bnd_particles_ops {
  psc_bnd_particles_ops_single() {
    name                    = "single";
    setup                   = psc_bnd_particles_sub<mparticles_t>::setup;
    unsetup                 = psc_bnd_particles_sub<mparticles_t>::unsetup;
    exchange_particles      = psc_bnd_particles_sub<mparticles_t>::exchange_particles;
    //open_calc_moments       = psc_bnd_particles_sub_open_calc_moments;
  }
} psc_bnd_particles_single_ops;
