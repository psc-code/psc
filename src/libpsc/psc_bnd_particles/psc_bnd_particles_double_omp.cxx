
#include "psc_bnd_particles_private.h"
#include "psc_particles_as_double.h"

#include "psc_bnd_particles_common_omp.cxx"

// ======================================================================
// psc_bnd_particles: subclass "double_omp"

struct psc_bnd_particles_ops psc_bnd_particles_double_omp_ops = {
  .name                    = "double_omp",
  .setup                   = psc_bnd_particles_sub_setup,
  .unsetup                 = psc_bnd_particles_sub_unsetup,
  .exchange_particles      = psc_bnd_particles_sub_exchange_particles,
};
