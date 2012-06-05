
#include "psc_push_fields_private.h"

static void
psc_push_fields_sub_push_E(struct psc_push_fields *push, struct psc_fields *flds_base)
{
}

static void
psc_push_fields_sub_push_H(struct psc_push_fields *push, struct psc_fields *flds_base)
{
}

// ======================================================================
// psc_push_fields: subclass "none"

struct psc_push_fields_ops psc_push_fields_none_ops = {
  .name                  = "none",
  .push_E                = psc_push_fields_sub_push_E,
  .push_H                = psc_push_fields_sub_push_H,
};
