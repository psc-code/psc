
#include "psc_bnd_fields_private.h"
#include "psc_fields_as_c.h"

#include "psc_bnd_fields_common.cxx"

// ======================================================================
// psc_bnd_fields: subclass "c"

struct psc_bnd_fields_ops psc_bnd_fields_c_ops = {
  .name                  = "c",
  .fill_ghosts_E         = psc_bnd_fields_sub_fill_ghosts_E,
  .fill_ghosts_H         = psc_bnd_fields_sub_fill_ghosts_H,
  .add_ghosts_J          = psc_bnd_fields_sub_add_ghosts_J,
};
