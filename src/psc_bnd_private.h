
#ifndef PSC_BND_PRIVATE_H
#define PSC_BND_PRIVATE_H

#include <psc_bnd.h>

struct psc_bnd {
  struct mrc_obj obj;
};

struct psc_bnd_ops {
  MRC_OBJ_OPS;
  void (*add_ghosts)(struct psc_bnd *bnd, mfields_base_t *flds, int mb, int me);
  void (*fill_ghosts)(struct psc_bnd *bnd, mfields_base_t *flds, int mb, int me);
  void (*exchange_particles)(struct psc_bnd *bnd, mparticles_base_t *particles);
};

// ======================================================================

extern struct psc_bnd_ops psc_bnd_c_ops;

#define to_psc_bnd(o) (container_of(o, struct psc_bnd, obj))
#define psc_bnd_ops(bnd) ((struct psc_bnd_ops *)((bnd)->obj.ops))

#endif
