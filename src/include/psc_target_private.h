
#ifndef PSC_TARGET_PRIVATE_H
#define PSC_TARGET_PRIVATE_H

#include <psc_target.h>

struct psc_target {
  struct mrc_obj obj;
};

struct psc_target_ops {
  MRC_SUBCLASS_OPS(struct psc_target);
  bool (*is_inside)(struct psc_target *target, double x[3]);
  void (*init_npt)(struct psc_target *target, int pop, double x[3],
		   struct psc_particle_npt *npt);
};

#define psc_target_ops(target) ((struct psc_target_ops *)((target)->obj.ops))

#endif


