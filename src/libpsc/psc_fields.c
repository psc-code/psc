
#include "psc.h"
#include "psc_fields_private.h"

#include <mrc_profile.h>

// ======================================================================
// psc_fields_init

static void
psc_fields_init()
{
  mrc_class_register_subclass(&mrc_class_psc_fields, &psc_fields_c_ops);
  mrc_class_register_subclass(&mrc_class_psc_fields, &psc_fields_single_ops);
  mrc_class_register_subclass(&mrc_class_psc_fields, &psc_fields_fortran_ops);
}

// ======================================================================
// psc_fields class

struct mrc_class_psc_fields mrc_class_psc_fields = {
  .name             = "psc_fields",
  .size             = sizeof(struct psc_fields),
  .init             = psc_fields_init,
};

