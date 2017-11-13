
#ifndef PSC_PARTICLES_SINGLE_BY_KIND_H
#define PSC_PARTICLES_SINGLE_BY_KIND_H

#include "psc_particles_private.h"

#include <psc.h>

#include "psc_particle_single_by_kind.h"

struct psc_mparticles_single_by_kind {
  struct bk_mparticles *bkmprts;
};

particle_single_by_kind_t *psc_mparticles_single_by_kind_get_one(struct psc_mparticles *mprts,
								 int p, unsigned int n);


#endif
