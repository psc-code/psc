
#ifndef PSC_PUSH_PARTICLES_SINGLE_H
#define PSC_PUSH_PARTICLES_SINGLE_H

#include "psc.h"

#include "psc_particles_as_single.h"
#include "psc_fields_as_c.h"

typedef fields_single_t fields_cache_t;
typedef fields_cache_t fields_curr_t;

void psc_push_particles_single_1vb_push_yz(struct psc_push_particles *push,
					   mparticles_base_t *particles_base,
					   mfields_base_t *flds_base);

#endif
