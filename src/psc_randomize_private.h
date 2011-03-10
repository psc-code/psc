
#ifndef PSC_RANDOMIZE_PRIVATE_H
#define PSC_RANDOMIZE_PRIVATE_H

#include <psc_randomize.h>

struct psc_randomize {
  struct mrc_obj obj;
};

struct psc_randomize_ops {
  MRC_OBJ_OPS;
  void (*run)(struct psc_randomize *randomize, mparticles_base_t *particles);
};

// ======================================================================

extern struct psc_randomize_ops psc_randomize_none_ops;
extern struct psc_randomize_ops psc_randomize_fortran_ops;

#define to_psc_randomize(o) (container_of(o, struct psc_randomize, obj))
#define psc_randomize_ops(randomize) ((struct psc_randomize_ops *)((randomize)->obj.ops))

#endif
