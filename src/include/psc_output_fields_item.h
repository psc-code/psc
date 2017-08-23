
#ifndef PSC_OUTPUT_FIELDS_ITEM_H
#define PSC_OUTPUT_FIELDS_ITEM_H

#include <mrc_obj.h>

#include "psc.h"

MRC_CLASS_DECLARE(psc_output_fields_item, struct psc_output_fields_item);

void psc_output_fields_item_set_psc_bnd(struct psc_output_fields_item *item,
					struct psc_bnd *bnd);
mfields_c_t *psc_output_fields_item_create_mfields(struct psc_output_fields_item *item);
void psc_output_fields_item_run(struct psc_output_fields_item *item,
				struct psc_mfields *flds, struct psc_mparticles *particles,
				mfields_c_t *res);

#endif

