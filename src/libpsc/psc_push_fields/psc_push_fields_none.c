
#include "psc_push_fields_private.h"

static void
psc_push_fields_sub_push_mflds_E(struct psc_push_fields *push, struct psc_mfields *mflds_base)
{
}

static void
psc_push_fields_sub_push_mflds_H(struct psc_push_fields *push, struct psc_mfields *mflds_base)
{
}

// ======================================================================
// psc_push_fields: subclass "none"

struct psc_push_fields_ops psc_push_fields_none_ops = {
  .name                  = "none",
  .push_mflds_E          = psc_push_fields_sub_push_mflds_E,
  .push_mflds_H          = psc_push_fields_sub_push_mflds_H,
};
