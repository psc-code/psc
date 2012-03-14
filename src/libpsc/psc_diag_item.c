
#include "psc_diag_item_private.h"

int
psc_diag_item_nr_values(struct psc_diag_item *item)
{
  return item->n_values;
}

const char *
psc_diag_item_title(struct psc_diag_item *item, int i)
{
  return item->names[i];
}

void
psc_diag_item_run(struct psc_diag_item *item, struct psc *psc,
		  double *result)
{
  item->run(psc, result);
}
