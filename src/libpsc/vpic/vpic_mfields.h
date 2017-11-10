
#ifndef VPIC_MFIELDS_H
#define VPIC_MFIELDS_H

#include "vpic_iface.h"

#include <vpic.h>

// ======================================================================
// vpic_mfields

struct vpic_mfields {
  field_array_t *field_array;

  vpic_mfields() :
    field_array(0) {
  }

  void clear_jf();
  void synchronize_jf();
};

#endif

