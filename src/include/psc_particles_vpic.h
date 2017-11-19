
#ifndef PSC_PARTICLES_VPIC_H
#define PSC_PARTICLES_VPIC_H

#include "psc_particles_private.h"
#include "psc_particles_single.h"

#include "../libpsc/vpic/vpic_iface.h" // FIXME path

#define PTYPE PTYPE_VPIC
#include "psc_particle_common.h"
#include "psc_particles_common.h"
#undef PTYPE

#endif
