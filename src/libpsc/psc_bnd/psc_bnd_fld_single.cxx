
#include "psc_bnd_fld.h"
#include "psc_bnd_private.h"
#include "psc_fields_as_single.h"

#define psc_bnd_fld_sub_copy_to_buf   psc_bnd_fld_single_copy_to_buf
#define psc_bnd_fld_sub_copy_from_buf psc_bnd_fld_single_copy_from_buf
#define psc_bnd_fld_sub_add_from_buf  psc_bnd_fld_single_add_from_buf

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

