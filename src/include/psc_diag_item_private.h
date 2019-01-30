
#ifndef PSC_DIAG_ITEM_PRIVATE_H
#define PSC_DIAG_ITEM_PRIVATE_H

#include "psc_diag_item.h"
#include "particles.hxx"
#include "fields3d.hxx"

struct psc_diag_item {
  struct mrc_obj obj;
};

struct psc_diag_item_ops {
  MRC_SUBCLASS_OPS(struct psc_diag_item);
  
  void (*run)(struct psc_diag_item *item,
	      MparticlesBase& mprts, MfieldsStateBase& mflds, double *result);
  int nr_values;
  const char *title[6]; // FIXME ugly hardcoded 6
};

#define psc_diag_item_ops(item) ((struct psc_diag_item_ops *)((item)->obj.ops))

#endif
