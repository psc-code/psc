
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
};

// ======================================================================

#define psc_push_fields_ops(push_fields) ((struct psc_push_fields_ops *)((push_fields)->obj.ops))

#endif
