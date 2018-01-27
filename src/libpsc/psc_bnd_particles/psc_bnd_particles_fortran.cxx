
#include "psc_bnd_particles_private.h"
#include "psc_particles_as_fortran.h"

#define NO_OPEN_BC
#include "psc_bnd_particles_common.cxx"

// ======================================================================
// psc_bnd_particles: subclass "fortran"

struct psc_bnd_particles_ops_fortran : psc_bnd_particles_ops {
  psc_bnd_particles_ops_fortran() {
    name                    = "fortran";
    setup                   = psc_bnd_particles_sub<mparticles_t>::setup;
    unsetup                 = psc_bnd_particles_sub<mparticles_t>::unsetup;
    exchange_particles      = psc_bnd_particles_sub<mparticles_t>::exchange_particles;
  }
} psc_bnd_particles_fortran_ops;
