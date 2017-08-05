
#ifndef PSC_PARTICLES_PRIVATE_H
#define PSC_PARTICLES_PRIVATE_H

#include "psc_particles.h"

// FIXME, hack to aid getting rid of ::n_part member
#define N_PART _n_part

struct psc_particles {
  struct mrc_obj obj;
  struct psc_mparticles *mprts;
  int N_PART;
  int n_alloced;
  int p; //< patch number
  unsigned int flags;
};

static inline int
psc_particles_size(struct psc_particles *prts)
{
  return prts->N_PART;
}

static inline void
psc_particles_resize(struct psc_particles *prts, int n_prts)
{
  assert(n_prts <= prts->n_alloced);
  prts->N_PART = n_prts;
}

static inline void
psc_particles_set_n_prts(struct psc_particles *prts, int n_prts)
{
  // FIXME, same as above w/o the assert, should go away...
  prts->N_PART = n_prts;
}

struct psc_particles_ops {
  MRC_SUBCLASS_OPS(struct psc_particles);
};

#define psc_particles_ops(prts) ((struct psc_particles_ops *) ((prts)->obj.ops))

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
