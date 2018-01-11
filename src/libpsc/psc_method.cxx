
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
// psc_method_setup_partition

void
psc_method_setup_partition(struct psc_method *method, struct psc *psc,
			   int *n_prts_by_patch)
{
  struct psc_method_ops *ops = psc_method_ops(method);
  assert(ops && ops->setup_partition);

  ops->setup_partition(method, psc, n_prts_by_patch);
}

// ----------------------------------------------------------------------
// psc_method_set_ic_particles

void
psc_method_set_ic_particles(struct psc_method *method, struct psc *psc,
			    int *n_prts_by_patch)
{
  struct psc_method_ops *ops = psc_method_ops(method);
  assert(ops && ops->set_ic_particles);

  ops->set_ic_particles(method, psc, n_prts_by_patch);
}

// ----------------------------------------------------------------------
// psc_method_set_ic_fields

void
psc_method_set_ic_fields(struct psc_method *method, struct psc *psc)
{
  struct psc_method_ops *ops = psc_method_ops(method);
  assert(ops && ops->set_ic_fields);

  ops->set_ic_fields(method, psc);
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

extern struct psc_method_ops psc_method_ops_default;
extern struct psc_method_ops psc_method_ops_vpic;

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

struct mrc_class_psc_method_ : mrc_class_psc_method {
  mrc_class_psc_method_() {
    name             = "psc_method";
    size             = sizeof(struct psc_method);
    init             = psc_method_init;
  }
} mrc_class_psc_method;




