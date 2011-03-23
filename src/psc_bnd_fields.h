
#ifndef PSC_BND_FIELDS_H
#define PSC_BND_FIELDS_H

#include <mrc_obj.h>

#include "psc.h"

MRC_CLASS_DECLARE(psc_bnd_fields, struct psc_bnd_fields);

void psc_bnd_fields_fill_ghosts_b_H(struct psc_bnd_fields *bnd, mfields_base_t *flds);

#endif
