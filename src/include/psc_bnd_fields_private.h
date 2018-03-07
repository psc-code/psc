
#ifndef PSC_BND_FIELDS_PRIVATE_H
#define PSC_BND_FIELDS_PRIVATE_H

#include <psc_bnd_fields.h>

struct psc_bnd_fields {
  struct mrc_obj obj;
};

struct psc_bnd_fields_ops {
  MRC_SUBCLASS_OPS(struct psc_bnd_fields);
};

// ======================================================================

#endif
