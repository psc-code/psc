
#include "psc_bnd_fld.h"
#include "psc_bnd_private.h"
#include "psc_particles_as_single.h"

// ======================================================================
// psc_bnd: subclass "single"

struct psc_bnd_ops psc_bnd_single_ops = {
  .name                    = "single",
  .create_ddc              = psc_bnd_fld_single_create,
  .add_ghosts              = psc_bnd_fld_single_add_ghosts,
  .fill_ghosts             = psc_bnd_fld_single_fill_ghosts,
};
