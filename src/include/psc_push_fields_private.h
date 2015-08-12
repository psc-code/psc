
#ifndef PSC_PUSH_FIELDS_PRIVATE_H
#define PSC_PUSH_FIELDS_PRIVATE_H

#include <psc_push_fields.h>

struct psc_push_fields {
  struct mrc_obj obj;
  // parameters
  int variant; //< 0: default, 1: optimized version with fewer fill_ghosts()

  // state
  struct psc_bnd_fields *bnd_fields;
};

struct psc_push_fields_ops {
  MRC_SUBCLASS_OPS(struct psc_push_fields);
  void (*push_E)(struct psc_push_fields *push, struct psc_fields *flds);
  void (*push_H)(struct psc_push_fields *push, struct psc_fields *flds);
  void (*push_mflds_E)(struct psc_push_fields *push, struct psc_mfields *mflds);
  void (*push_mflds_H)(struct psc_push_fields *push, struct psc_mfields *mflds);
  void (*pml_a)(struct psc_push_fields *push, struct psc_fields *flds);
  void (*pml_b)(struct psc_push_fields *push, struct psc_fields *flds);
};

// ======================================================================

extern struct psc_push_fields_ops psc_push_fields_auto_ops;
extern struct psc_push_fields_ops psc_push_fields_c_ops;
extern struct psc_push_fields_ops psc_push_fields_single_ops;
extern struct psc_push_fields_ops psc_push_fields_fortran_ops;
extern struct psc_push_fields_ops psc_push_fields_cbe_ops;
extern struct psc_push_fields_ops psc_push_fields_cuda_ops;
extern struct psc_push_fields_ops psc_push_fields_cuda2_ops;
extern struct psc_push_fields_ops psc_push_fields_acc_ops;
extern struct psc_push_fields_ops psc_push_fields_mix_ops;
extern struct psc_push_fields_ops psc_push_fields_none_ops;

#define psc_push_fields_ops(push_fields) ((struct psc_push_fields_ops *)((push_fields)->obj.ops))

#endif
