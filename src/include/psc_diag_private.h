
#ifndef PSC_DIAG_PRIVATE_H
#define PSC_DIAG_PRIVATE_H

#include <psc_diag.h>

#include <stdio.h>

struct psc_diag {
  struct mrc_obj obj;

  // parameters
  int every_step;

  // internal
  std::vector<psc_diag_item*> items_;
  FILE *file_;
  int rank_;
};

#endif
