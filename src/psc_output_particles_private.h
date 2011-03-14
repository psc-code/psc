
#ifndef PSC_OUTPUT_PARTICLES_PRIVATE_H
#define PSC_OUTPUT_PARTICLES_PRIVATE_H

#include <psc_output_particles.h>

struct psc_output_particles {
  struct mrc_obj obj;
};

struct psc_output_particles_ops {
  MRC_SUBCLASS_OPS(struct psc_output_particles);
  void (*run)(struct psc_output_particles *output_particles,
	      mparticles_base_t *particles);
};

// ======================================================================

extern struct psc_output_particles_ops psc_output_particles_fortran_ops;

#define psc_output_particles_ops(output_particles) ((struct psc_output_particles_ops *)((output_particles)->obj.ops))

#endif
