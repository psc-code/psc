
#include "psc_bnd_fields_private.h"
#include "psc_fields_c.h"
#include "psc_fields_single.h"

#include "psc_bnd_fields_common.cxx"

// ======================================================================
// psc_bnd_fields: subclass "c"

using bnd_fields_ops_c = bnd_fields_ops<mfields_c_t>;

struct psc_bnd_fields_ops psc_bnd_fields_c_ops = {
  .name                  = "c",
  .fill_ghosts_E         = bnd_fields_ops_c::fill_ghosts_E,
  .fill_ghosts_H         = bnd_fields_ops_c::fill_ghosts_H,
  .add_ghosts_J          = bnd_fields_ops_c::add_ghosts_J,
};

// ======================================================================
// psc_bnd_fields: subclass "single"

using bnd_fields_ops_single = bnd_fields_ops<mfields_single_t>;

struct psc_bnd_fields_ops psc_bnd_fields_single_ops = {
  .name                  = "single",
  .fill_ghosts_E         = bnd_fields_ops_single::fill_ghosts_E,
  .fill_ghosts_H         = bnd_fields_ops_single::fill_ghosts_H,
  .add_ghosts_J          = bnd_fields_ops_single::add_ghosts_J,
};
