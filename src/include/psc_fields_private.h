
#ifndef PSC_FIELDS_PRIVATE_H
#define PSC_FIELDS_PRIVATE_H

#include "psc_fields.h"

struct psc_fields {
  struct mrc_obj obj;
  void *data;
};

struct psc_fields_ops {
  MRC_SUBCLASS_OPS(struct psc_fields);
};

#define psc_fields_ops(prts) ((struct psc_fields_ops *) ((prts)->obj.ops))

// ======================================================================

extern struct psc_fields_ops psc_fields_c_ops;
extern struct psc_fields_ops psc_fields_single_ops;
extern struct psc_fields_ops psc_fields_fortran_ops;
extern struct psc_fields_ops psc_fields_cuda_ops;
extern struct psc_fields_ops psc_fields_cuda2_ops;
extern struct psc_fields_ops psc_fields_acc_ops;

extern struct psc_mfields_ops psc_mfields_cuda2_ops;
extern struct psc_mfields_ops psc_mfields_acc_ops;

#endif
