
#include "psc_bnd_particles_private.h"
#include "psc_particles_as_fortran.h"

#define NO_OPEN_BC
#include "psc_bnd_particles_common.cxx"

// ======================================================================
// psc_bnd_particles: subclass "fortran"

struct psc_bnd_particles_ops_fortran : psc_bnd_particles_ops {
  psc_bnd_particles_ops_fortran() {
    name                    = "fortran";
    setup                   = psc_bnd_particles_sub_setup;
    unsetup                 = psc_bnd_particles_sub_unsetup;
    exchange_mprts_prep     = mparticles_ddcp<mparticles_t>::exchange_mprts_prep;
    exchange_mprts_post     = mparticles_ddcp<mparticles_t>::exchange_mprts_post;
    exchange_particles      = psc_bnd_particles_sub_exchange_particles;
  }
} psc_bnd_particles_fortran_ops;
