
#include "psc_bnd_fields_private.h"

#include "psc.h"
#include "psc_glue.h"
#include <mrc_profile.h>

// ----------------------------------------------------------------------
// psc_bnd_fields_none_fill_ghosts_b_H

static void
psc_bnd_fields_none_fill_ghosts_b_H(struct psc_bnd_fields *bnd,
				    mfields_base_t *flds_base)
{
  // FIXME, should check that no pulses are set, either,
  // or better, move the pulses -> fortran subclass
  for (int d = 0; d < 3; d++) {
    assert(psc.domain.bnd_fld_lo[d] == BND_FLD_PERIODIC);
    assert(psc.domain.bnd_fld_hi[d] == BND_FLD_PERIODIC);
  }
}

// ======================================================================
// psc_bnd_fields: subclass "none"

struct psc_bnd_fields_ops psc_bnd_fields_none_ops = {
  .name                  = "none",
  .fill_ghosts_b_H       = psc_bnd_fields_none_fill_ghosts_b_H,
};
