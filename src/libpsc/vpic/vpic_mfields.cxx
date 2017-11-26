
#include "vpic_iface.h"

#include <cassert>

// ======================================================================
// vpic_mfields

// ----------------------------------------------------------------------
// C wrappers

void vpic_mfields_accumulate_rho_p(FieldArray *vmflds, Particles *vmprts)
{
  for (Particles::Iter sp = vmprts->begin(); sp != vmprts->end(); ++sp) {
    TIC accumulate_rho_p(vmflds, &*sp); TOC( accumulate_rho_p, 1);
  }
}

