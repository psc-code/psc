
#include "psc_bnd_fld.h"
#include "psc_bnd_private.h"

// ======================================================================
// psc_bnd: subclass "c"

struct psc_bnd_ops psc_bnd_c_ops = {
  .name                    = "c",
  .create_ddc              = psc_bnd_fld_c_create,
  .add_ghosts              = psc_bnd_fld_c_add_ghosts,
  .fill_ghosts             = psc_bnd_fld_c_fill_ghosts,
};
