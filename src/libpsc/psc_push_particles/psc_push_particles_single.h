
#ifndef PSC_PUSH_PARTICLES_SINGLE_H
#define PSC_PUSH_PARTICLES_SINGLE_H

#include "psc.h"

#include "psc_particles_as_single.h"
#include "psc_fields_as_c.h"

void psc_push_particles_single_1vb_push_a_yz(struct psc_push_particles *push,
					     struct psc_particles *particles_base,
					     struct psc_fields *flds_base);
void psc_push_particles_single2_1vb_push_a_yz(struct psc_push_particles *push,
					      struct psc_particles *particles_base,
					      struct psc_fields *flds_base);

#endif
