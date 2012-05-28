
#include "psc_bnd_c.h"
#include "psc_bnd_private.h"
#include "ddc_particles.h"
#include "psc_particles_as_single.h"

#include "psc_bnd_common2.c"

// ======================================================================
// psc_bnd: subclass "single2"

struct psc_bnd_ops psc_bnd_single2_ops = {
  .name                  = "single2",
  .size                  = sizeof(struct psc_bnd_sub),
  .setup                 = psc_bnd_sub_setup,
  .unsetup               = psc_bnd_sub_unsetup,
  .exchange_particles    = psc_bnd_sub_exchange_particles,

  .create_ddc            = psc_bnd_fields_c_create,
  .add_ghosts            = psc_bnd_fields_c_add_ghosts,
  .fill_ghosts           = psc_bnd_fields_c_fill_ghosts,
};
