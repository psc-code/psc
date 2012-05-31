#ifndef PSC_PUSH_FIELDS_H
#define PSC_PUSH_FIELDS_H


#include <mrc_obj.h>

#include "psc.h"

MRC_CLASS_DECLARE(psc_push_fields, struct psc_push_fields);

struct psc_bnd_fields *psc_push_fields_get_bnd_fields(struct psc_push_fields *push);

void psc_push_fields_step_a(struct psc_push_fields *push, mfields_base_t *flds);
void psc_push_fields_step_b_H(struct psc_push_fields *push, mfields_base_t *flds);
void psc_push_fields_step_b_E(struct psc_push_fields *push, mfields_base_t *flds);

#endif
