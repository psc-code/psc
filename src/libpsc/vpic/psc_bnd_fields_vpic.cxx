
#include "psc_bnd_fields_private.h"

#include "vpic_iface.h"

// ----------------------------------------------------------------------
// psc_bnd_fields: subclass "vpic"

struct psc_bnd_fields_ops_vpic : psc_bnd_fields_ops {
  psc_bnd_fields_ops_vpic() {
    name                  = "vpic";
  }
} psc_bnd_fields_vpic_ops;

