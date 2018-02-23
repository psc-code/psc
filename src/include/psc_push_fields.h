#ifndef PSC_PUSH_FIELDS_H
#define PSC_PUSH_FIELDS_H

#include <mrc_obj.h>

#include "psc.h"

BEGIN_C_DECLS

MRC_CLASS_DECLARE(psc_push_fields, struct psc_push_fields);

struct psc_bnd_fields *psc_push_fields_get_bnd_fields(struct psc_push_fields *push);

void psc_push_fields_push_H(struct psc_push_fields *push, struct psc_mfields *flds,
			    double dt_fac);
void psc_push_fields_push_E(struct psc_push_fields *push, struct psc_mfields *flds,
			    double dt_fac);

END_C_DECLS


#endif
