
#ifndef VPIC_MPARTICLES_H
#define VPIC_MPARTICLES_H

#include "vpic_iface.h"

#include "simulation.h"
#include <vpic.h>

// ======================================================================
// vpic_mparticles

struct vpic_mparticles {
  vpic_mparticles(species_t*& sl) : p_(sl) { }
  Particles p_;
};

#endif

