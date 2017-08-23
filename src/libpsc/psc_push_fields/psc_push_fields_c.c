
#include "psc_push_fields_private.h"
#include "psc_fields_as_c.h"

#include "psc_push_fields_common.c"

// ======================================================================
// psc_push_fields: subclass "c"

struct psc_push_fields_ops psc_push_fields_c_ops = {
  .name                  = "c",
  .push_mflds_E          = psc_push_fields_sub_push_mflds_E,
  .push_mflds_H          = psc_push_fields_sub_push_mflds_H,
};
