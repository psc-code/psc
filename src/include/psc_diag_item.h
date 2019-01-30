
#ifndef PSC_DIAG_ITEM_H
#define PSC_DIAG_ITEM_H

#include <mrc_obj.h>

#include "psc.h"
#include "particles.hxx"
#include "fields3d.hxx"

MRC_CLASS_DECLARE(psc_diag_item, struct psc_diag_item);

int  psc_diag_item_nr_values(struct psc_diag_item *item);
const char *psc_diag_item_title(struct psc_diag_item *item, int i);
void psc_diag_item_run(struct psc_diag_item *item,
		       MparticlesBase& mprts, MfieldsStateBase& mflds,
		       double *result);
		       

#endif
