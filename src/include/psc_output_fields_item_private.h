
#ifndef PSC_OUTPUT_FIELDS_ITEM_PRIVATE_H
#define PSC_OUTPUT_FIELDS_ITEM_PRIVATE_H

#include <psc_output_fields_item.h>

struct psc_output_fields_item {
  struct mrc_obj obj;
  struct psc_bnd *bnd;
};

enum {
  POFI_ADD_GHOSTS = 1, // this item needs to have ghost points added to interior points
  POFI_BY_KIND    = 2, // this item needs to be replicated by kind
};

#define POFI_MAX_COMPS (16)

struct psc_output_fields_item_ops {
  MRC_SUBCLASS_OPS(struct psc_output_fields_item);
  void (*run_all)(struct psc_output_fields_item *item,
		  struct psc_mfields *mflds, struct psc_mparticles *mprts,
		  struct psc_mfields *mres);
  int nr_comp;
  const char *fld_names[POFI_MAX_COMPS];
  unsigned int flags;
};

#define psc_output_fields_item_ops(item)			\
  ((struct psc_output_fields_item_ops *)((item)->obj.ops))

// ======================================================================

#endif
