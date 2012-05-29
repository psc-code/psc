
#ifndef PSC_BND_PRIVATE_H
#define PSC_BND_PRIVATE_H

#include <psc_push_fields.h>

struct psc_push_fields {
  struct mrc_obj obj;
  struct psc_bnd_fields *bnd_fields;
};

struct psc_push_fields_ops {
  MRC_SUBCLASS_OPS(struct psc_push_fields);
  void (*push_a_E)(struct psc_push_fields *push, struct psc_fields *flds);
  void (*push_a_H)(struct psc_push_fields *push, struct psc_fields *flds);
  void (*push_b_H)(struct psc_push_fields *push, struct psc_fields *flds);
  void (*push_b_E)(struct psc_push_fields *push, struct psc_fields *flds);
  void (*pml_a)(struct psc_push_fields *push, struct psc_fields *flds);
  void (*pml_b)(struct psc_push_fields *push, struct psc_fields *flds);
};

// ======================================================================

extern struct psc_push_fields_ops psc_push_fields_c_ops;
extern struct psc_push_fields_ops psc_push_fields_single_ops;
extern struct psc_push_fields_ops psc_push_fields_fortran_ops;
extern struct psc_push_fields_ops psc_push_fields_cbe_ops;
extern struct psc_push_fields_ops psc_push_fields_cuda_ops;

#define psc_push_fields_ops(push_fields) ((struct psc_push_fields_ops *)((push_fields)->obj.ops))

#endif
