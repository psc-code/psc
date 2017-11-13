
#include "psc_method_private.h"

// ======================================================================
// psc_method

// ----------------------------------------------------------------------
// psc_method_do_setup

void
psc_method_do_setup(struct psc_method *method, struct psc *psc)
{
  struct psc_method_ops *ops = psc_method_ops(method);
  assert(ops && ops->do_setup);

  ops->do_setup(method, psc);
}

// ----------------------------------------------------------------------
// psc_method_setup_partition_and_particles

void
psc_method_setup_partition_and_particles(struct psc_method *method, struct psc *psc)
{
  struct psc_method_ops *ops = psc_method_ops(method);
  assert(ops && ops->setup_partition_and_particles);

  ops->setup_partition_and_particles(method, psc);
}

// ----------------------------------------------------------------------
// psc_method_setup_fields

void
psc_method_setup_fields(struct psc_method *method, struct psc *psc)
{
  struct psc_method_ops *ops = psc_method_ops(method);
  assert(ops && ops->setup_fields);

  ops->setup_fields(method, psc);
}

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




