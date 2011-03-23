
#include "psc_push_fields_private.h"
#include <mrc_profile.h>

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
    static int pr_e, pr_h;
    if (!pr_e) {
      pr_e = prof_register("push_fields_a_e", 1., 0, 0);
      pr_h = prof_register("push_fields_a_h", 1., 0, 0);
    }
    assert(ops->push_a_e);
    assert(ops->push_a_h);

    prof_start(pr_e);
    ops->push_a_e(push, flds);
    prof_stop(pr_e);

    prof_start(pr_h);
    ops->push_a_h(push, flds);
    prof_stop(pr_h);
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
    static int pr_e, pr_h;
    if (!pr_e) {
      pr_e = prof_register("push_fields_b_e", 1., 0, 0);
      pr_h = prof_register("push_fields_b_h", 1., 0, 0);
    }
    assert(ops->push_b_h);
    assert(ops->push_b_e);

    prof_start(pr_e);
    ops->push_b_h(push, flds);
    prof_stop(pr_e);

    prof_start(pr_h);
    ops->push_b_e(push, flds);
    prof_stop(pr_h);
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

