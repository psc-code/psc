
#ifndef PSC_METHOD_PRIVATE_H
#define PSC_METHOD_PRIVATE_H

#include "psc_method.h"

struct psc_method {
  struct mrc_obj obj;
};

struct psc_method_ops {
  MRC_SUBCLASS_OPS(struct psc_method);
};

#define psc_method_ops(method) ((struct psc_method_ops *)((method)->obj.ops))

#endif
