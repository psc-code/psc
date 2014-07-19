#ifndef PSC_MARDER_H
#define PSC_MARDER_H

#include <mrc_obj.h>

#include "psc.h"

MRC_CLASS_DECLARE(psc_marder, struct psc_marder);

void psc_marder_run(struct psc_marder *marder, struct psc_mfields *mflds,
		    struct psc_mparticles *mprts);

#endif
