
#ifndef PSC_DIAG_PRIVATE_H
#define PSC_DIAG_PRIVATE_H

#include <psc_diag.h>

#include <stdio.h>

struct psc_diag {
  struct mrc_obj obj;
  int every_step;
  struct psc_diag_item **items;

  FILE *file;
};

#endif
