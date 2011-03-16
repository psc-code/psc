
#include "psc_push_fields_private.h"

// ======================================================================
// forward to subclass

void
psc_push_fields_step_a(struct psc_push_fields *push, mfields_base_t *flds)
{
  struct psc_push_fields_ops *ops = psc_push_fields_ops(push);
  assert(ops->step_a);
  ops->step_a(push, flds);
}

void
psc_push_fields_step_b(struct psc_push_fields *push, mfields_base_t *flds)
{
  struct psc_push_fields_ops *ops = psc_push_fields_ops(push);
  assert(ops->step_b);
  ops->step_b(push, flds);
}

// ======================================================================
// psc_push_fields_init

static void
psc_push_fields_init()
{
  mrc_class_register_subclass(&mrc_class_psc_push_fields, &psc_push_fields_c_ops);
  mrc_class_register_subclass(&mrc_class_psc_push_fields, &psc_push_fields_fortran_ops);
#ifdef USE_CBE
  mrc_class_register_subclass(&mrc_class_psc_push_fields, &psc_push_fields_cbe_ops);
#endif
}

// ======================================================================
// psc_push_fields class

struct mrc_class_psc_push_fields mrc_class_psc_push_fields = {
  .name             = "psc_push_fields",
  .size             = sizeof(struct psc_push_fields),
  .init             = psc_push_fields_init,
};

