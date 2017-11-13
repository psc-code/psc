
#ifndef PSC_PARTICLE_SINGLE_BY_KIND_H
#define PSC_PARTICLE_SINGLE_BY_KIND_H

#include "psc_particles_private.h"

#include <psc.h>

#define PTYPE PTYPE_SINGLE_BY_KIND
#include "psc_particle_buf_common.h"
#undef PTYPE

struct psc_mparticles_single_by_kind {
  struct psc_mparticles_single_by_kind_patch *patch;

  struct bk_mparticles *bkmprts;
};

particle_single_by_kind_t *psc_mparticles_single_by_kind_get_one(struct psc_mparticles *mprts,
								 int p, unsigned int n);

void psc_mparticles_single_by_kind_patch_reserve(struct psc_mparticles *mprts, int p,
						 unsigned int new_capacity);


#endif
