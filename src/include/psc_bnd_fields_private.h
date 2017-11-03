
#ifndef PSC_BND_FIELDS_PRIVATE_H
#define PSC_BND_FIELDS_PRIVATE_H

#include <psc_bnd_fields.h>

struct psc_bnd_fields {
  struct mrc_obj obj;
};

struct psc_bnd_fields_ops {
  MRC_SUBCLASS_OPS(struct psc_bnd_fields);
  void (*fill_ghosts_E)(struct psc_bnd_fields *bnd, struct psc_mfields *mflds);
  void (*fill_ghosts_H)(struct psc_bnd_fields *bnd, struct psc_mfields *mflds);
  void (*add_ghosts_J)(struct psc_bnd_fields *bnd, struct psc_mfields *mflds);
};

// ======================================================================

extern struct psc_bnd_fields_ops psc_bnd_fields_fortran_ops;
extern struct psc_bnd_fields_ops psc_bnd_fields_none_ops;
extern struct psc_bnd_fields_ops psc_bnd_fields_c_ops;
extern struct psc_bnd_fields_ops psc_bnd_fields_single_ops;
extern struct psc_bnd_fields_ops psc_bnd_fields_cuda_ops;

#define psc_bnd_fields_ops(bnd_fields) ((struct psc_bnd_fields_ops *)((bnd_fields)->obj.ops))

#endif
