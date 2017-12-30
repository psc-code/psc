
#ifndef PSC_CHECKS_H
#define PSC_CHECKS_H

#include <mrc_obj.h>

#include "psc.h"

BEGIN_C_DECLS

MRC_CLASS_DECLARE(psc_checks, struct psc_checks);

void psc_checks_gauss(struct psc_checks *checks, struct psc *psc);

void psc_checks_continuity_before_particle_push(struct psc_checks *checks, struct psc *psc);
void psc_checks_continuity_after_particle_push(struct psc_checks *checks, struct psc *psc);

END_C_DECLS

#endif
