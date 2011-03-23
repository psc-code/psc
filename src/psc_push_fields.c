
#include "psc_push_fields_private.h"

// ======================================================================
// forward to subclass

void
psc_push_fields_step_a(struct psc_push_fields *push, mfields_base_t *flds)
{
  struct psc_push_fields_ops *ops = psc_push_fields_ops(push);
  if (psc.domain.use_pml) {
    // FIXME, pml routines sehould be split into E, H push + ghost points, too
    // pml could become a separate push_fields subclass
    assert(ops->pml_a);
    ops->pml_a(push, flds);
  } else {
    assert(ops->step_a);
    ops->step_a(push, flds);
  }
}

void
psc_push_fields_step_b(struct psc_push_fields *push, mfields_base_t *flds)
{
  struct psc_push_fields_ops *ops = psc_push_fields_ops(push);
  if (psc.domain.use_pml) {
    assert(ops->pml_b);
    ops->pml_b(push, flds);
  } else {
    assert(ops->step_b);
    ops->step_b(push, flds);
  }
}

// ======================================================================
// psc_push_fields_init

static void
psc_push_fields_init()
{
  mrc_class_register_subclass(&mrc_class_psc_push_fields, &psc_push_fields_c_ops);
  mrc_class_register_subclass(&mrc_class_psc_push_fields, &psc_push_fields_fortran_ops);
}

// ======================================================================
// psc_push_fields class

struct mrc_class_psc_push_fields mrc_class_psc_push_fields = {
  .name             = "psc_push_fields",
  .size             = sizeof(struct psc_push_fields),
  .init             = psc_push_fields_init,
};

