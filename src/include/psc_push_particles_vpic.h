
#ifndef PSC_PUSH_PARTICLES_VPIC_H
#define PSC_PUSH_PARTICLES_VPIC_H

#include "psc_push_particles_private.h"

// ======================================================================
// psc_push_particles_vpic

struct psc_push_particles_vpic {
  struct vpic_push_particles *vpushp;
};

#define psc_push_particles_vpic(push) mrc_to_subobj(push, struct psc_push_particles_vpic)

#endif
