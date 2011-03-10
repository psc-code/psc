
#include <mrc_obj.h>

#include "psc.h"

extern struct mrc_class mrc_class_psc_output_fields;

MRC_OBJ_DEFINE_STANDARD_METHODS(psc_output_fields, struct psc_output_fields)

void psc_output_fields_run(struct psc_output_fields *output_fields,
			   mfields_base_t *flds, mparticles_base_t *particles);
