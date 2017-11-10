
#ifndef PSC_COLLISION_PRIVATE_H
#define PSC_COLLISION_PRIVATE_H

#include <psc_collision.h>

struct psc_collision {
  struct mrc_obj obj;
};

struct psc_collision_ops {
  MRC_SUBCLASS_OPS(struct psc_collision);
  void (*run)(struct psc_collision *collision, struct psc_mparticles *mprts);
};

// ======================================================================

extern struct psc_collision_ops psc_collision_none_ops;
extern struct psc_collision_ops psc_collision_fortran_ops;
extern struct psc_collision_ops psc_collision_single_ops;
extern struct psc_collision_ops psc_collision_vpic_ops;

#define psc_collision_ops(collision) ((struct psc_collision_ops *)((collision)->obj.ops))

#endif
