
#ifndef PSC_HEATING_H
#define PSC_HEATING_H

#include "psc.h"

MRC_CLASS_DECLARE(psc_heating, struct psc_heating);

void psc_heating_run(struct psc_heating *heating, struct psc_mparticles *mprts_base,
		     struct psc_mfields *mflds_base);

#endif
