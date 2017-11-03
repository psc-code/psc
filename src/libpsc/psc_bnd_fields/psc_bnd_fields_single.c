
#include "psc_bnd_fields_private.h"
#include "psc_fields_as_single.h"

#include "psc_bnd_fields_common.c"

// ======================================================================
// psc_bnd_fields: subclass "single"

struct psc_bnd_fields_ops psc_bnd_fields_single_ops = {
  .name                  = "single",
  .fill_ghosts_E         = psc_bnd_fields_sub_fill_ghosts_E,
  .fill_ghosts_H         = psc_bnd_fields_sub_fill_ghosts_H,
  .add_ghosts_J          = psc_bnd_fields_sub_add_ghosts_J,
};
