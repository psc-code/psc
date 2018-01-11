
#ifndef PSC_SORT_PRIVATE_H
#define PSC_SORT_PRIVATE_H

#include <psc_sort.h>

struct psc_sort {
  struct mrc_obj obj;
  int every; //< sort every so many steps
};

struct psc_sort_ops {
  MRC_SUBCLASS_OPS(struct psc_sort);
  void (*run)(struct psc_sort *sort, struct psc_mparticles *mprts);
};

// ======================================================================

#define psc_sort_ops(sort) ((struct psc_sort_ops *)((sort)->obj.ops))

#endif
