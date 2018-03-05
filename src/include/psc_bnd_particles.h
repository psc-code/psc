
#ifndef PSC_BND_PARTICLES_H
#define PSC_BND_PARTICLES_H

#include <mrc_obj.h>

#include "psc.h"

MRC_CLASS_DECLARE(psc_bnd_particles, struct psc_bnd_particles);

void psc_bnd_particles_set_psc(struct psc_bnd_particles *bnd, struct psc *psc);

#endif
