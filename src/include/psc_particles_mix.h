
#ifndef PSC_PARTICLE_MIX_H
#define PSC_PARTICLE_MIX_H

#include "psc_particles_private.h"

struct psc_mparticles_mix {
  struct psc_mparticles *sub;
};

#define psc_mparticles_mix(mprts) mrc_to_subobj(mprts, struct psc_mparticles_mix)

#endif
