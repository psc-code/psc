
#include "psc_bnd_particles_private.h"
#include "../psc_bnd/ddc_particles.h"
#include "psc_particles_as_double.h"

#include "psc_bnd_particles_common.c"

// ======================================================================
// psc_bnd_particles: subclass "double"

struct psc_bnd_particles_ops psc_bnd_particles_double_ops = {
  .name                    = "double",
  .setup                   = psc_bnd_particles_sub_setup,
  .unsetup                 = psc_bnd_particles_sub_unsetup,
  .exchange_particles      = psc_bnd_particles_sub_exchange_particles,
  .exchange_particles_prep = psc_bnd_particles_sub_exchange_particles_prep,
  .exchange_particles_post = psc_bnd_particles_sub_exchange_particles_post,
  .open_calc_moments       = psc_bnd_particles_sub_open_calc_moments,
};
