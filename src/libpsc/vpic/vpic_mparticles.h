
#ifndef VPIC_MPARTICLES_H
#define VPIC_MPARTICLES_H

#include "vpic_iface.h"

#include "simulation.h"
#include <vpic.h>

// ======================================================================
// vpic_mparticles

struct vpic_mparticles : Particles {
  vpic_mparticles(species_t*& sl) : Particles(sl) { }
};

#endif

