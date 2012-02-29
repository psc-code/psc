
#ifndef PSC_DIAG_H
#define PSC_DIAG_H

#include <mrc_obj.h>

#include "psc.h"

struct psc_diag_item {
  void (*run)(struct psc *psc, double *result);
  int n_values;
  const char *names[];
};

extern struct psc_diag_item psc_diag_item_em_energy;
extern struct psc_diag_item psc_diag_item_particle_energy;

MRC_CLASS_DECLARE(psc_diag, struct psc_diag);

void psc_diag_set_items(struct psc_diag *diag, struct psc_diag_item **items);
void psc_diag_run(struct psc_diag *diag, struct psc *psc);

#endif
