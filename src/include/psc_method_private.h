
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

BEGIN_C_DECLS

// maybe useful for non-default subclasses

void psc_method_default_output(struct psc_method *method, struct psc *psc,
			       int stats_every,
			       MfieldsBase& mflds, MparticlesBase& mprts);

END_C_DECLS

#endif
