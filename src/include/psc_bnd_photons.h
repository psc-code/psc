
#ifndef PSC_BND_PHOTONS_H
#define PSC_BND_PHOTONS_H

#include <mrc_obj.h>

#include "psc.h"

MRC_CLASS_DECLARE(psc_bnd_photons, struct psc_bnd_photons);

void psc_bnd_photons_set_psc(struct psc_bnd_photons *bnd, struct psc *psc);
void psc_bnd_photons_exchange(struct psc_bnd_photons *bnd, mphotons_t *mphotons);

#endif
