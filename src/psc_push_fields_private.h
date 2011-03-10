
#ifndef PSC_BND_PRIVATE_H
#define PSC_BND_PRIVATE_H

#include <psc_push_fields.h>

struct psc_push_fields {
  struct mrc_obj obj;
};

struct psc_push_fields_ops {
  MRC_OBJ_OPS;
  void (*step_a)(struct psc_push_fields *push, mfields_base_t *flds);
  void (*step_b)(struct psc_push_fields *push, mfields_base_t *flds);
};

// ======================================================================

extern struct psc_push_fields_ops psc_push_fields_c_ops;
extern struct psc_push_fields_ops psc_push_fields_fortran_ops;

#define to_psc_push_fields(o) (container_of(o, struct psc_push_fields, obj))
#define psc_push_fields_ops(push_fields) ((struct psc_push_fields_ops *)((push_fields)->obj.ops))

#endif
