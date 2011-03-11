
#ifndef PSC_SORT_PRIVATE_H
#define PSC_SORT_PRIVATE_H

#include <psc_sort.h>

struct psc_sort {
  struct mrc_obj obj;
};

struct psc_sort_ops {
  MRC_SUBCLASS_OPS(struct psc_sort);
  void (*run)(struct psc_sort *sort, mparticles_base_t *particles);
};

// ======================================================================

extern struct psc_sort_ops psc_sort_none_ops;
extern struct psc_sort_ops psc_sort_fortran_ops;
extern struct psc_sort_ops psc_sort_qsort_ops;
extern struct psc_sort_ops psc_sort_countsort_ops;
extern struct psc_sort_ops psc_sort_countsort2_ops;

#define to_psc_sort(o) (container_of(o, struct psc_sort, obj))
#define psc_sort_ops(sort) ((struct psc_sort_ops *)((sort)->obj.ops))

#endif
