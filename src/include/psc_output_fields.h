
#include <mrc_obj.h>

#include "psc.h"

MRC_CLASS_DECLARE(psc_output_fields, struct psc_output_fields);

BEGIN_C_DECLS

void psc_output_fields_set_psc(struct psc_output_fields *output_fields,
			       struct psc *psc);
void psc_output_fields_run(struct psc_output_fields *output_fields,
			   struct psc_mfields *flds, struct psc_mparticles *particles);

END_C_DECLS
