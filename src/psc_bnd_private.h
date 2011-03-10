
#ifndef PSC_BND_PRIVATE_H
#define PSC_BND_PRIVATE_H

#include <psc_bnd.h>

struct psc_bnd {
  struct mrc_obj obj;
};

struct psc_bnd_ops {
  MRC_OBJ_OPS;
};

// ======================================================================

struct psc_bnd_c {
  struct mrc_ddc *ddc;
  struct ddc_particles *ddcp;
};

extern struct psc_bnd_ops psc_bnd_c_ops;


#endif
