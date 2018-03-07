
#ifndef PSC_BND_PRIVATE_H
#define PSC_BND_PRIVATE_H

#include <psc_bnd.h>

struct psc_bnd {
  struct mrc_obj obj;
  struct psc *psc;
  struct mrc_ddc *ddc;
};

struct psc_bnd_ops {
  MRC_SUBCLASS_OPS(struct psc_bnd);

  // for field exchange
  void (*create_ddc)(struct psc_bnd *bnd);
  void (*add_ghosts)(struct psc_bnd *bnd, struct psc_mfields *flds, int mb, int me);
  void (*fill_ghosts)(struct psc_bnd *bnd, struct psc_mfields *flds, int mb, int me);
  void (*reset)(struct psc_bnd *bnd);
};

// ======================================================================

#define psc_bnd_ops(bnd) ((struct psc_bnd_ops *)((bnd)->obj.ops))

#endif
