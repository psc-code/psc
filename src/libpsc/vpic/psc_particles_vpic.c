
#include "psc_particles_vpic.h"

#include "vpic_iface.h"

// ----------------------------------------------------------------------
// psc_mparticles: subclass "vpic"
  
struct psc_mparticles_ops psc_mparticles_vpic_ops = {
  .name                    = "vpic",
  .size                    = sizeof(struct psc_mparticles_vpic),
};

