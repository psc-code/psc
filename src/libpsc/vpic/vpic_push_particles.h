
#ifndef VPIC_PUSH_PARTICLES_H
#define VPIC_PUSH_PARTICLES_H

#include "vpic_iface.h"

#include <vpic.h>

// ======================================================================
// vpic_push_particles

struct vpic_push_particles {
  interpolator_array_t *interpolator_array;
  accumulator_array_t *accumulator_array;
};

void vpic_clear_accumulator_array(struct vpic_push_particles *vpushp);
void vpic_advance_p(struct vpic_push_particles *vpushp, struct vpic_mparticles *vmprts);
void vpic_reduce_accumulator_array(struct vpic_push_particles *vpushp);
void vpic_boundary_p(struct vpic_push_particles *vpushp, struct vpic_mparticles *vmprts,
		     struct vpic_mfields *vmflds);
void vpic_push_particles_unload_accumulator_array(struct vpic_push_particles *vpushp,
						  struct vpic_mfields *vmflds);
void vpic_push_particles_load_interpolator_array(struct vpic_push_particles *vpushp,
						 struct vpic_mfields *vmflds);

#endif

