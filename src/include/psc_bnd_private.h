
#ifndef PSC_BND_PRIVATE_H
#define PSC_BND_PRIVATE_H

#include <psc_bnd.h>

struct psc_bnd {
  struct mrc_obj obj;
  struct psc *psc;
};

struct psc_bnd_ops {
  MRC_SUBCLASS_OPS(struct psc_bnd);
};

// ======================================================================

#define psc_bnd_ops(bnd) ((struct psc_bnd_ops *)((bnd)->obj.ops))

#endif
