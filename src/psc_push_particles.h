
#include <mrc_obj.h>

#include "psc.h"

extern struct mrc_class mrc_class_psc_push_particles;

MRC_OBJ_DEFINE_STANDARD_METHODS(psc_push_particles, struct psc_push_particles)

void psc_push_particles_run(struct psc_push_particles *push,
			    mparticles_base_t *particles, mfields_base_t *flds);
void psc_push_particles_push_yz_a(struct psc_push_particles *push,
				  mparticles_base_t *particles, mfields_base_t *flds);
void psc_push_particles_push_yz_b(struct psc_push_particles *push,
				  mparticles_base_t *particles, mfields_base_t *flds);
