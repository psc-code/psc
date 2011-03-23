
#ifndef PSC_BND_FIELDS_PRIVATE_H
#define PSC_BND_FIELDS_PRIVATE_H

#include <psc_bnd_fields.h>

struct psc_bnd_fields {
  struct mrc_obj obj;
};

struct psc_bnd_fields_ops {
  MRC_SUBCLASS_OPS(struct psc_bnd_fields);
  void (*fill_ghosts_b_H)(struct psc_bnd_fields *bnd, mfields_base_t *flds);
};

// ======================================================================

extern struct psc_bnd_fields_ops psc_bnd_fields_fortran_ops;

#define psc_bnd_fields_ops(bnd_fields) ((struct psc_bnd_fields_ops *)((bnd_fields)->obj.ops))

#endif
