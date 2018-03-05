
#include "psc_bnd_particles_private.h"

#include "vpic_iface.h"

// ----------------------------------------------------------------------
// psc_bnd_particles: subclass "vpic"

struct psc_bnd_particles_ops_vpic : psc_bnd_particles_ops {
  psc_bnd_particles_ops_vpic() {
    name                  = "vpic";
  }
} psc_bnd_particles_vpic_ops;

