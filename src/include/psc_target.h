
#ifndef PSC_TARGET_H
#define PSC_TARGET_H

#include <psc.h>

BEGIN_C_DECLS

MRC_CLASS_DECLARE(psc_target, struct psc_target);

bool psc_target_is_inside(struct psc_target *target, double x[3]);
void psc_target_init_npt(struct psc_target *target, int pop, double x[3],
			 struct psc_particle_npt *npt);

END_C_DECLS

#endif

