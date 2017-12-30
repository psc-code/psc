
#ifndef PSC_BND_H
#define PSC_BND_H

#include <mrc_obj.h>

#include "psc.h"

BEGIN_C_DECLS

MRC_CLASS_DECLARE(psc_bnd, struct psc_bnd);

void psc_bnd_set_psc(struct psc_bnd *bnd, struct psc *psc);
void psc_bnd_add_ghosts(struct psc_bnd *bnd, struct psc_mfields *flds, int mb, int me);
void psc_bnd_fill_ghosts(struct psc_bnd *bnd, struct psc_mfields *flds, int mb, int me);

END_C_DECLS

#endif
