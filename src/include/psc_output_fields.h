
#include <mrc_obj.h>

#include "psc.h"

MRC_CLASS_DECLARE(psc_output_fields, struct psc_output_fields);

void psc_output_fields_run(struct psc_output_fields *output_fields,
			   mfields_base_t *flds, mparticles_base_t *particles);
