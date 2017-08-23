#ifndef PSC_PUSH_FIELDS_H
#define PSC_PUSH_FIELDS_H

#include <mrc_obj.h>

#include "psc.h"

MRC_CLASS_DECLARE(psc_push_fields, struct psc_push_fields);

struct psc_bnd_fields *psc_push_fields_get_bnd_fields(struct psc_push_fields *push);

void psc_push_fields_step_a(struct psc_push_fields *push, struct psc_mfields *flds);
void psc_push_fields_step_b1(struct psc_push_fields *push, struct psc_mfields *flds);
void psc_push_fields_step_b2(struct psc_push_fields *push, struct psc_mfields *flds);

#endif
