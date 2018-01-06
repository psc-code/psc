
#ifndef PSC_INJECT_H
#define PSC_INJECT_H

#include <psc.h>

BEGIN_C_DECLS

MRC_CLASS_DECLARE(psc_inject, struct psc_inject);

void psc_inject_run(struct psc_inject *inject, struct psc_mparticles *mprts_base,
		    struct psc_mfields *mflds_base);

END_C_DECLS

#endif

