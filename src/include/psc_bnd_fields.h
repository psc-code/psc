
#ifndef PSC_BND_FIELDS_H
#define PSC_BND_FIELDS_H

#include <mrc_obj.h>

#include "psc.h"

MRC_CLASS_DECLARE(psc_bnd_fields, struct psc_bnd_fields);

void psc_bnd_fields_fill_ghosts_E(struct psc_bnd_fields *bnd, struct psc_mfields *mflds);
void psc_bnd_fields_fill_ghosts_H(struct psc_bnd_fields *bnd, struct psc_mfields *mflds);
void psc_bnd_fields_add_ghosts_J(struct psc_bnd_fields *bnd, struct psc_mfields *mflds);

#endif
