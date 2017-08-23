
#include <mrc_obj.h>

#include "psc.h"

MRC_CLASS_DECLARE(psc_output_fields_collection, struct psc_output_fields_collection);

void psc_output_fields_collection_set_psc(struct psc_output_fields_collection *output_fields_collection,
					  struct psc *psc);
void psc_output_fields_collection_run(struct psc_output_fields_collection *output_fields_collection,
				      struct psc_mfields *flds, struct psc_mparticles *particles);
