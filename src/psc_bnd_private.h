
#ifndef PSC_BND_PRIVATE_H
#define PSC_BND_PRIVATE_H

#include <psc_bnd.h>

struct psc_bnd {
  struct mrc_obj obj;
};

struct psc_bnd_ops {
  MRC_OBJ_OPS;
};

void psc_bnd_register(struct psc_bnd_ops *ops);

void psc_bnd_register_c(void);

#endif
