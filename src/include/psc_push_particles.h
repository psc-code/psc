
#ifndef PSC_PUSH_PARTICLES_H
#define PSC_PUSH_PARTICLES_H

#include <mrc_obj.h>

#include "psc.h"

BEGIN_C_DECLS

MRC_CLASS_DECLARE(psc_push_particles, struct psc_push_particles);

void psc_push_particles_prep(struct psc_push_particles *push,
			     struct psc_mparticles *mprts, struct psc_mfields *mflds);
void psc_push_particles_run(struct psc_push_particles *push,
			    struct psc_mparticles *mprts, struct psc_mfields *mflds);
void psc_push_particles_stagger(struct psc_push_particles *push,
				struct psc_mparticles *mprts, struct psc_mfields *mflds);
unsigned int psc_push_particles_get_mp_flags(struct psc_push_particles *push);

END_C_DECLS

#endif
