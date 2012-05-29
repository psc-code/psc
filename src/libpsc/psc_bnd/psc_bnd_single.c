
#include "psc_bnd_fld.h"
#include "psc_bnd_private.h"
#include "ddc_particles.h"
#include "psc_particles_as_single.h"

#include "psc_bnd_common.c"

// ======================================================================
// psc_bnd: subclass "single"

struct psc_bnd_ops psc_bnd_single_ops = {
  .name                    = "single",
  .size                    = sizeof(struct psc_bnd_sub),
  .setup                   = psc_bnd_sub_setup,
  .unsetup                 = psc_bnd_sub_unsetup,
  .exchange_particles      = psc_bnd_sub_exchange_particles,
  .exchange_particles_prep = psc_bnd_sub_exchange_particles_prep,
  .exchange_particles_post = psc_bnd_sub_exchange_particles_post,

  .create_ddc              = psc_bnd_fld_single_create,
  .add_ghosts              = psc_bnd_fld_single_add_ghosts,
  .fill_ghosts             = psc_bnd_fld_single_fill_ghosts,
};
