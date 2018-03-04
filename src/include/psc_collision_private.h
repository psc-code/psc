
#ifndef PSC_COLLISION_PRIVATE_H
#define PSC_COLLISION_PRIVATE_H

#include <psc_collision.h>

struct psc_collision {
  struct mrc_obj obj;
  int every;
  double nu;
};

struct psc_collision_ops {
  MRC_SUBCLASS_OPS(struct psc_collision);
};

// ======================================================================

#define psc_collision_ops(collision) ((struct psc_collision_ops *)((collision)->obj.ops))

#endif
