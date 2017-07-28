
#ifndef PSC_PARTICLES_PRIVATE_H
#define PSC_PARTICLES_PRIVATE_H

#include "psc_particles.h"

// FIXME, hack to aid getting rid of ::n_part member
#define N_PART n_part

struct psc_particles {
  struct mrc_obj obj;
  struct psc_mparticles *mprts;
  int N_PART;
  int p; //< patch number
  unsigned int flags;
};

struct psc_particles_ops {
  MRC_SUBCLASS_OPS(struct psc_particles);
  void (*reorder)(struct psc_particles *prts);
};

#define psc_particles_ops(prts) ((struct psc_particles_ops *) ((prts)->obj.ops))

typedef void (*psc_particles_copy_to_func_t)(struct psc_particles *,
					     struct psc_particles *,
					     unsigned int);
typedef void (*psc_particles_copy_from_func_t)(struct psc_particles *,
					       struct psc_particles *,
					       unsigned int);

// ======================================================================

extern struct psc_particles_ops psc_particles_c_ops;
extern struct psc_particles_ops psc_particles_single_ops;
extern struct psc_particles_ops psc_particles_double_ops;
extern struct psc_particles_ops psc_particles_single_by_block_ops;
extern struct psc_particles_ops psc_particles_fortran_ops;
extern struct psc_particles_ops psc_particles_cuda_ops;
extern struct psc_particles_ops psc_particles_cuda2_ops;
extern struct psc_particles_ops psc_particles_acc_ops;

extern struct psc_mparticles_ops psc_mparticles_cuda2_ops;
extern struct psc_mparticles_ops psc_mparticles_acc_ops;

#endif
