
#ifndef PSC_PUSH_PARTICLES_DOUBLE_H
#define PSC_PUSH_PARTICLES_DOUBLE_H

#include "psc.h"

#include "psc_particles_as_double.h"
#include "psc_fields_as_c.h"

void psc_push_particles_1vb_double_push_a_yz(struct psc_push_particles *push,
					     struct psc_particles *prts_base,
					     struct psc_fields *flds_base);
void psc_push_particles_1vbec_double_push_a_yz(struct psc_push_particles *push,
					       struct psc_particles *prts_base,
					       struct psc_fields *flds_base);

#endif
