
#ifndef PSC_OUTPUT_FIELDS_ITEM_PRIVATE_H
#define PSC_OUTPUT_FIELDS_ITEM_PRIVATE_H

#include <psc_output_fields_item.h>

struct psc_output_fields_item {
  struct mrc_obj obj;
  char *name;
  int nr_comp;
  char *fld_names[6];
  void (*calc)(struct psc *psc, mfields_base_t *flds, mparticles_base_t *particles,
	       mfields_c_t *f);
};

// ======================================================================

#endif
