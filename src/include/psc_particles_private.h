
#ifndef PSC_PARTICLES_PRIVATE_H
#define PSC_PARTICLES_PRIVATE_H

#include "psc_particles.h"

struct psc_particles {
  struct mrc_obj obj;
  struct psc_mparticles *mprts;
  int p; //< patch number
};

#endif
