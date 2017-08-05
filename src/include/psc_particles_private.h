
#ifndef PSC_PARTICLES_PRIVATE_H
#define PSC_PARTICLES_PRIVATE_H

#include "psc_particles.h"

struct psc_particles {
  struct mrc_obj obj;
  struct psc_mparticles *mprts;
  int p; //< patch number
};

static inline int
psc_particles_size(struct psc_particles *prts)
{
  return psc_mparticles_n_prts_by_patch(prts->mprts, prts->p);
}

static inline void
psc_particles_resize(struct psc_particles *prts, int n_prts)
{
  psc_mparticles_resize_patch(prts->mprts, prts->p, n_prts);
}

static inline void
psc_particles_set_n_prts(struct psc_particles *prts, int n_prts)
{
  psc_mparticles_set_n_prts_by_patch(prts->mprts, prts->p, n_prts);
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
