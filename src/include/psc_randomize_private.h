
#ifndef PSC_RANDOMIZE_PRIVATE_H
#define PSC_RANDOMIZE_PRIVATE_H

#include <psc_randomize.h>

struct psc_randomize {
  struct mrc_obj obj;
};

struct psc_randomize_ops {
  MRC_SUBCLASS_OPS(struct psc_randomize);
  void (*run)(struct psc_randomize *randomize, struct psc_mparticles *mprts_base);
};

// ======================================================================

extern struct psc_randomize_ops psc_randomize_none_ops;
extern struct psc_randomize_ops psc_randomize_c_ops;
extern struct psc_randomize_ops psc_randomize_fortran_ops;

#define psc_randomize_ops(randomize) ((struct psc_randomize_ops *)((randomize)->obj.ops))

#endif
