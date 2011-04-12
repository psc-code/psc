
#include <mrc_obj.h>

#include "psc.h"

MRC_CLASS_DECLARE(psc_output_particles, struct psc_output_particles);

void psc_output_particles_run(struct psc_output_particles *output_particles,
			      mparticles_base_t *particles);
