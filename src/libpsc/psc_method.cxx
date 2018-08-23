
#include "psc_method_private.h"

// ======================================================================
// psc_method

// ----------------------------------------------------------------------
// psc_method_init

extern struct psc_method_ops psc_method_ops_default;

static void
psc_method_init(void)
{
  mrc_class_register_subclass(&mrc_class_psc_method, &psc_method_ops_default);
}

// ----------------------------------------------------------------------
// psc_method class

struct mrc_class_psc_method_ : mrc_class_psc_method {
  mrc_class_psc_method_() {
    name             = "psc_method";
    size             = sizeof(struct psc_method);
    init             = psc_method_init;
  }
} mrc_class_psc_method;




