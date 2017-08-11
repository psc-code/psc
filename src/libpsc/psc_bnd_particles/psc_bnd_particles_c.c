
#include "psc_bnd_particles_private.h"
#include "psc_particles_as_c.h"

#define NO_OPEN_BC
#include "psc_bnd_particles_common.c"

// ======================================================================
// psc_bnd_particles: subclass "c"

struct psc_bnd_particles_ops psc_bnd_particles_c_ops = {
  .name                    = "c",
  .setup                   = psc_bnd_particles_sub_setup,
  .unsetup                 = psc_bnd_particles_sub_unsetup,
  .exchange_particles      = psc_bnd_particles_sub_exchange_particles,
  .open_calc_moments       = psc_bnd_particles_sub_open_calc_moments,
};
