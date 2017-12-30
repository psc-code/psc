
#ifndef PSC_BND_PARTICLES_H
#define PSC_BND_PARTICLES_H

#include <mrc_obj.h>

#include "psc.h"

BEGIN_C_DECLS

MRC_CLASS_DECLARE(psc_bnd_particles, struct psc_bnd_particles);

void psc_bnd_particles_set_psc(struct psc_bnd_particles *bnd, struct psc *psc);
void psc_bnd_particles_exchange(struct psc_bnd_particles *bnd, struct psc_mparticles *mprts);

END_C_DECLS

#endif
