
#ifndef PSC_BND_PHOTONS_PRIVATE_H
#define PSC_BND_PHOTONS_PRIVATE_H

#include <psc_bnd_photons.h>

struct psc_bnd_photons {
  struct mrc_obj obj;
  struct psc *psc;
  struct ddc_particles *ddcp;
};

struct psc_bnd_photons_ops {
  MRC_SUBCLASS_OPS(struct psc_bnd);
};

#endif
