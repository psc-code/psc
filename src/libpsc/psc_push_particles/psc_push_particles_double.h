
#ifndef PSC_PUSH_PARTICLES_DOUBLE_H
#define PSC_PUSH_PARTICLES_DOUBLE_H

#include "psc.h"

#include "psc_particles_as_double.h"
#include "psc_fields_as_c.h"

typedef struct psc_fields fields_cache_t;
typedef struct psc_fields fields_curr_t;

void psc_push_particles_double_1vb_push_yz(struct psc_push_particles *push,
					   mparticles_base_t *particles_base,
					   mfields_base_t *flds_base);

#endif
