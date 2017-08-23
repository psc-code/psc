
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
  void (*add_ghosts_prep)(struct psc_bnd *bnd, struct psc_fields *pf, int mb, int me);
  void (*add_ghosts_post)(struct psc_bnd *bnd, struct psc_fields *pf, int mb, int me);
  void (*fill_ghosts)(struct psc_bnd *bnd, struct psc_mfields *flds, int mb, int me);
  void (*fill_ghosts_prep)(struct psc_bnd *bnd, struct psc_fields *pf, int mb, int me);
  void (*fill_ghosts_post)(struct psc_bnd *bnd, struct psc_fields *pf, int mb, int me);
};

// ======================================================================

extern struct psc_bnd_ops psc_bnd_auto_ops;
extern struct psc_bnd_ops psc_bnd_c_ops;
extern struct psc_bnd_ops psc_bnd_single_ops;
extern struct psc_bnd_ops psc_bnd_mix_ops;
extern struct psc_bnd_ops psc_bnd_cuda_ops;

#define psc_bnd_ops(bnd) ((struct psc_bnd_ops *)((bnd)->obj.ops))

#endif
