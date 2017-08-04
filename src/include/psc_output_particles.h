
#include <mrc_obj.h>

#include "psc.h"

#ifndef PSC_OUTPUT_PARTICLES_H
#define PSC_OUTPUT_PARTICLES_H

MRC_CLASS_DECLARE(psc_output_particles, struct psc_output_particles);

void psc_output_particles_run(struct psc_output_particles *output_particles,
			      struct psc_mparticles *particles);
#endif			      
