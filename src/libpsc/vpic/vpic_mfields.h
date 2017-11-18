
#ifndef VPIC_MFIELDS_H
#define VPIC_MFIELDS_H

#include "vpic_iface.h"

#include "field_array.h"
#include "hydro_array.h"

#include <vpic.h>

// ======================================================================
// vpic_mfields

struct vpic_mfields : FieldArray {
  vpic_mfields(Grid g, MaterialList m_list, float damp)
    : FieldArray(g, m_list, damp) { }
  
  void accumulate_rho_p(vpic_mparticles *vmprts);
};

struct vpic_mfields_hydro : HydroArray {
};

#endif

