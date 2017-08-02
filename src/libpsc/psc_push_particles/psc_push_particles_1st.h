
#ifndef PSC_PUSH_PARTICLES_1ST_H
#define PSC_PUSH_PARTICLES_1ST_H

#include "psc.h"
#include "psc_fields_as_c.h"
#include "psc_particles_as_c.h"

void psc_push_particles_1st_push_mprts_xz(struct psc_push_particles *push,
					  struct psc_mparticles *mprts,
					  struct psc_mfields *mflds);
void psc_push_particles_1st_push_mprts_yz(struct psc_push_particles *push,
					  struct psc_mparticles *mprts,
					  struct psc_mfields *mflds);

void psc_push_particles_1vb_c_push_mprts_yz(struct psc_push_particles *push,
					    struct psc_mparticles *mprts,
					    struct psc_mfields *mflds);

void psc_push_particles_1sff_push_mprts_xz(struct psc_push_particles *push,
					   struct psc_mparticles *mprts,
					   struct psc_mfields *mflds);

#endif
