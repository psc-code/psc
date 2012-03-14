
#ifndef PSC_DIAG_PRIVATE_H
#define PSC_DIAG_PRIVATE_H

#include <psc_diag.h>

#include <stdio.h>

// FIXME, make the items children instead
#define MAX_ITEMS 20 // FIXME

struct psc_diag {
  struct mrc_obj obj;
  // parameters
  const char *items;
  int every_step;

  // internal
  int nr_items;
  struct psc_diag_item *item_list[MAX_ITEMS];
  FILE *file;
};

#endif
