
#ifndef PSC_PARTICLES_PRIVATE_H
#define PSC_PARTICLES_PRIVATE_H

#include "psc_particles.h"

struct psc_particles {
  struct mrc_obj obj;
  int n_part;
  unsigned int flags;
};

struct psc_particles_ops {
  MRC_SUBCLASS_OPS(struct psc_particles);
};

#define psc_particles_ops(prts) ((struct psc_particles_ops *) ((prts)->obj.ops))

// ======================================================================

extern struct psc_particles_ops psc_particles_single_ops;
extern struct psc_particles_ops psc_particles_double_ops;
extern struct psc_particles_ops psc_particles_fortran_ops;

#endif
