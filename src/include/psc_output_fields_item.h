
#ifndef PSC_OUTPUT_FIELDS_ITEM_H
#define PSC_OUTPUT_FIELDS_ITEM_H

#include <mrc_obj.h>

#include "psc.h"

MRC_CLASS_DECLARE(psc_output_fields_item, struct psc_output_fields_item);

void psc_output_fields_item_run(struct psc_output_fields_item *item,
				mfields_base_t *flds, mparticles_base_t *particles,
				mfields_c_t *res);

#endif

