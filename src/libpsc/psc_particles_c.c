
#include "psc.h"
#include "psc_particles_as_c.h"

#include "psc_particles_common.c"

// ======================================================================
// psc_mparticles: subclass "c"
  
struct psc_mparticles_ops psc_mparticles_c_ops = {
  .name                    = "c",
};

