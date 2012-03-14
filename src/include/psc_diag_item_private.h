
#ifndef PSC_DIAG_ITEM_PRIVATE_H
#define PSC_DIAG_ITEM_PRIVATE_H

#include "psc.h"

struct psc_diag_item {
  void (*run)(struct psc *psc, double *result);
  int n_values;
  const char *names[];
};

extern struct psc_diag_item psc_diag_item_em_energy;
extern struct psc_diag_item psc_diag_item_particle_energy;

#endif
