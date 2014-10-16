
#include "psc_bnd_particles_private.h"
#include "../psc_bnd/ddc_particles.h"
#include "psc_particles_as_fortran.h"

#include "psc_bnd_particles_common.c"

// ======================================================================
// psc_bnd_particles: subclass "fortran"

struct psc_bnd_particles_ops psc_bnd_particles_fortran_ops = {
  .name                    = "fortran",
  .setup                   = psc_bnd_particles_sub_setup,
  .unsetup                 = psc_bnd_particles_sub_unsetup,
  .exchange_particles      = psc_bnd_particles_sub_exchange_particles,
  .exchange_particles_prep = psc_bnd_particles_sub_exchange_particles_prep,
  .exchange_particles_post = psc_bnd_particles_sub_exchange_particles_post,
};
