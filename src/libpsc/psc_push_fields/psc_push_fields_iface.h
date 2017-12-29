
#ifndef PSC_PUSH_FIELDS_IFACE_H
#define PSC_PUSH_FIELDS_IFACE_H

#include <mrc_common.h>

#include <psc_fields_single.h>

BEGIN_C_DECLS

void psc_push_fields_single_push_E_xz(struct psc_push_fields *push, fields_single_t flds,
				      struct psc *psc, double dt_fac);
void psc_push_fields_single_push_H_xz(struct psc_push_fields *push, fields_single_t flds,
				      struct psc *psc, double dt_fac);

END_C_DECLS

#endif
