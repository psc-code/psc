
#ifndef VPIC_MFIELDS_H
#define VPIC_MFIELDS_H

#include "vpic_iface.h"

#include <vpic.h>

// ======================================================================
// vpic_mfields

struct vpic_mfields {
  field_array_t *field_array;
};

void vpic_mfields_clear_jf(struct vpic_mfields *vmflds);
void vpic_mfields_synchronize_jf(struct vpic_mfields *vmflds);

#endif

