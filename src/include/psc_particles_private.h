
#ifndef PSC_PARTICLES_PRIVATE_H
#define PSC_PARTICLES_PRIVATE_H

#include "psc_particles.h"

#include "psc_particles_double.h"

struct psc_particles {
  struct mrc_obj obj;
  int n_part;
  int n_alloced;
  unsigned int flags;
  struct psc_particle_double *particles;
};

struct psc_particles_ops {
  MRC_SUBCLASS_OPS(struct psc_particles);
};

// ======================================================================
// type "double"

static inline particle_double_t *
particles_double_get_one(struct psc_particles *prts, int n)
{
  return &prts->particles[n];
}

// ======================================================================

extern struct psc_particles_ops psc_particles_double_ops;

#endif
