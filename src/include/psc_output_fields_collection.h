
#ifndef PSC_OUTPUT_FIELDS_COLLECTION_H
#define PSC_OUTPUT_FIELDS_COLLECTION_H

#include <mrc_obj.h>

#include "psc.h"

MRC_CLASS_DECLARE(psc_output_fields_collection, struct psc_output_fields_collection);

BEGIN_C_DECLS

void psc_output_fields_collection_set_psc(struct psc_output_fields_collection *output_fields_collection,
					  struct psc *psc);
void psc_output_fields_collection_run(struct psc_output_fields_collection *output_fields_collection,
				      struct psc_mfields *flds, struct psc_mparticles *particles);

END_C_DECLS

#endif
