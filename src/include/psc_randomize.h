
#ifndef PSC_RANDOMIZE_H
#define PSC_RANDOMIZE_H

#include <mrc_obj.h>

#include "psc.h"

BEGIN_C_DECLS

MRC_CLASS_DECLARE(psc_randomize, struct psc_randomize);

void psc_randomize_run(struct psc_randomize *randomize, struct psc_mparticles *particles);

END_C_DECLS

#endif
