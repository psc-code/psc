
#ifndef PSC_DIAG_H
#define PSC_DIAG_H

#include <mrc_obj.h>

#include "psc.h"
#include "psc_diag_item.h"

MRC_CLASS_DECLARE(psc_diag, struct psc_diag);

void psc_diag_run(struct psc_diag *diag, MparticlesBase& mprts, MfieldsStateBase& mflds);

#endif
