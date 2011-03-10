
#include <mrc_obj.h>

#include "psc.h"

extern struct mrc_class mrc_class_psc_push_fields;

MRC_OBJ_DEFINE_STANDARD_METHODS(psc_push_fields, struct psc_push_fields)

void psc_push_fields_step_a(struct psc_push_fields *push, mfields_base_t *flds);
void psc_push_fields_step_b(struct psc_push_fields *push, mfields_base_t *flds);
