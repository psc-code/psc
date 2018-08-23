
#include "psc_method_private.h"

#include <setup_particles.hxx>
#include <psc_particles_double.h>

#include <vpic_iface.h>

// ======================================================================
// psc_method "vpic"

struct psc_method_vpic {
};

// ----------------------------------------------------------------------
// psc_method "vpic"

struct psc_method_ops_vpic : psc_method_ops {
  psc_method_ops_vpic() {
    name                          = "vpic";
    size                          = sizeof(struct psc_method_vpic);
  }
} psc_method_ops_vpic;
