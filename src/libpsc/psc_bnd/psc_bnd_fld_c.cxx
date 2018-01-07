
#include "psc_bnd_fld.h"
#include "psc_bnd_private.h"
#include "psc_fields_as_c.h"

#define psc_bnd_fld_sub_copy_to_buf   psc_bnd_fld_c_copy_to_buf
#define psc_bnd_fld_sub_copy_from_buf psc_bnd_fld_c_copy_from_buf
#define psc_bnd_fld_sub_add_from_buf  psc_bnd_fld_c_add_from_buf

#include "psc_bnd_fld.cxx"

// ======================================================================
// psc_bnd: subclass "c"

using psc_bnd_fld_ops_c = psc_bnd_fld_ops<mfields_c_t>;

struct psc_bnd_ops psc_bnd_c_ops = {
  .name                    = "c",
  .create_ddc              = psc_bnd_fld_ops_c::create,
  .add_ghosts              = psc_bnd_fld_ops_c::add_ghosts,
  .fill_ghosts             = psc_bnd_fld_ops_c::fill_ghosts,
};


