
#include "psc_bnd_fld.h"
#include "psc_bnd_private.h"
#include "psc_fields_as_single.h"

#include "psc_bnd_fld.cxx"

// ======================================================================
// psc_bnd: subclass "single"

using psc_bnd_fld_ops_single = psc_bnd_fld_ops<mfields_single_t>;

struct psc_bnd_ops psc_bnd_single_ops = {
  .name                    = "single",
  .create_ddc              = psc_bnd_fld_ops_single::create,
  .add_ghosts              = psc_bnd_fld_ops_single::add_ghosts,
  .fill_ghosts             = psc_bnd_fld_ops_single::fill_ghosts,
};

