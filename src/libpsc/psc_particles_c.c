
#include "psc.h"
#include "psc_particles_as_c.h"

// ======================================================================
// psc_mparticles: subclass "c"
  
#include "psc_particles_common.c"

struct psc_mparticles_ops psc_mparticles_c_ops = {
  .name                    = "c",
};

