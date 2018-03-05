
#ifndef PSC_BND_H
#define PSC_BND_H

#include <mrc_obj.h>

#include "psc.h"

MRC_CLASS_DECLARE(psc_bnd, struct psc_bnd);

void psc_bnd_set_psc(struct psc_bnd *bnd, struct psc *psc);

void psc_bnd_check_domain(struct psc_bnd *bnd); // FIXME

#endif
