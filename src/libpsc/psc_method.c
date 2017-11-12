
#include "psc_method_private.h"

// ======================================================================
// psc_method

// ----------------------------------------------------------------------
// psc_method_initialize

void
psc_method_initialize(struct psc_method *method, struct psc *psc)
{
  struct psc_method_ops *ops = psc_method_ops(method);
  assert(ops && ops->initialize);

  ops->initialize(method, psc);
}

// ----------------------------------------------------------------------
// psc_method_output

void
psc_method_output(struct psc_method *method, struct psc *psc)
{
  struct psc_method_ops *ops = psc_method_ops(method);
  assert(ops && ops->output);

  ops->output(method, psc);
}

// ----------------------------------------------------------------------
// psc_method_init

static void
psc_method_init(void)
{
  mrc_class_register_subclass(&mrc_class_psc_method, &psc_method_ops_default);
#ifdef USE_VPIC
  mrc_class_register_subclass(&mrc_class_psc_method, &psc_method_ops_vpic);
#endif
}

// ----------------------------------------------------------------------
// psc_method class

struct mrc_class_psc_method mrc_class_psc_method = {
  .name             = "psc_method",
  .size             = sizeof(struct psc_method),
  .init             = psc_method_init,
};




