
#ifndef PSC_INJECT_H
#define PSC_INJECT_H

#include <psc.h>

MRC_CLASS_DECLARE(psc_inject, struct psc_inject);

void psc_inject_run(struct psc_inject *inject, struct psc_mparticles *mprts_base,
		    struct psc_mfields *mflds_base);
#endif

