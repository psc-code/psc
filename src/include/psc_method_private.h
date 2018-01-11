
#ifndef PSC_METHOD_PRIVATE_H
#define PSC_METHOD_PRIVATE_H

#include "psc_method.h"

struct psc_method {
  struct mrc_obj obj;
};

struct psc_method_ops {
  MRC_SUBCLASS_OPS(struct psc_method);
  void (*do_setup)(struct psc_method *method, struct psc *psc);
  void (*setup_partition)(struct psc_method *method, struct psc *psc, int *n_prts_by_patch);
  void (*set_ic_particles)(struct psc_method *method, struct psc *psc, int *n_prts_by_patch);
  void (*set_ic_fields)(struct psc_method *method, struct psc *psc);
  void (*initialize)(struct psc_method *method, struct psc *psc);
  void (*output)(struct psc_method *method, struct psc *psc);
};

#define psc_method_ops(method) ((struct psc_method_ops *)((method)->obj.ops))

BEGIN_C_DECLS

// maybe useful for non-default subclasses

void psc_method_default_output(struct psc_method *method, struct psc *psc);

END_C_DECLS

#endif
