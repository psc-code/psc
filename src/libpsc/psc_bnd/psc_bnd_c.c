
#include "psc_bnd_c.h"
#include "psc_bnd_private.h"
#include "ddc_particles.h"
#include "psc_particles_as_c.h"

#include "psc_bnd_common.c"

// ======================================================================
// psc_bnd: subclass "c"

struct psc_bnd_ops psc_bnd_c_ops = {
  .name                  = "c",
  .size                  = sizeof(struct psc_bnd_sub),
  .setup                 = psc_bnd_sub_setup,
  .unsetup               = psc_bnd_sub_unsetup,
  .destroy               = psc_bnd_sub_destroy,
  .add_ghosts            = psc_bnd_sub_add_ghosts,
  .fill_ghosts           = psc_bnd_sub_fill_ghosts,
  .exchange_particles    = psc_bnd_sub_exchange_particles,
  .exchange_photons      = psc_bnd_sub_exchange_photons,
};
