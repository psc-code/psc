
#include "psc_push_fields_private.h"
#include <mrc_profile.h>

// ======================================================================
// forward to subclass

static inline void
psc_push_fields_push_a_e(struct psc_push_fields *push, mfields_base_t *flds)
{
  struct psc_push_fields_ops *ops = psc_push_fields_ops(push);
  static int pr;
  if (!pr) {
    pr = prof_register("push_fields_a_e", 1., 0, 0);
  }
  assert(ops->push_a_e);
  
  prof_start(pr);
  ops->push_a_e(push, flds);
  prof_stop(pr);
}

static inline void
psc_push_fields_push_a_h(struct psc_push_fields *push, mfields_base_t *flds)
{
  struct psc_push_fields_ops *ops = psc_push_fields_ops(push);
  static int pr;
  if (!pr) {
    pr = prof_register("push_fields_a_h", 1., 0, 0);
  }
  assert(ops->push_a_h);
  
  prof_start(pr);
  ops->push_a_h(push, flds);
  prof_stop(pr);
}

static inline void
psc_push_fields_push_b_e(struct psc_push_fields *push, mfields_base_t *flds)
{
  struct psc_push_fields_ops *ops = psc_push_fields_ops(push);
  static int pr;
  if (!pr) {
    pr = prof_register("push_fields_b_e", 1., 0, 0);
  }
  assert(ops->push_b_e);
  
  prof_start(pr);
  ops->push_b_e(push, flds);
  prof_stop(pr);
}

static inline void
psc_push_fields_push_b_h(struct psc_push_fields *push, mfields_base_t *flds)
{
  struct psc_push_fields_ops *ops = psc_push_fields_ops(push);
  static int pr;
  if (!pr) {
    pr = prof_register("push_fields_b_h", 1., 0, 0);
  }
  assert(ops->push_b_h);
  
  prof_start(pr);
  ops->push_b_h(push, flds);
  prof_stop(pr);
}

void
psc_push_fields_step_a(struct psc_push_fields *push, mfields_base_t *flds)
{
  if (psc.domain.use_pml) {
    // FIXME, pml routines sehould be split into E, H push + ghost points, too
    // pml could become a separate push_fields subclass
    struct psc_push_fields_ops *ops = psc_push_fields_ops(push);
    assert(ops->pml_a);
    ops->pml_a(push, flds);
  } else {
    psc_push_fields_push_a_e(push, flds);
    psc_push_fields_push_a_h(push, flds);
  }
}

void
psc_push_fields_step_b(struct psc_push_fields *push, mfields_base_t *flds)
{
  if (psc.domain.use_pml) {
    struct psc_push_fields_ops *ops = psc_push_fields_ops(push);
    assert(ops->pml_b);
    ops->pml_b(push, flds);
  } else {
    psc_push_fields_push_b_h(push, flds);
    psc_push_fields_push_b_e(push, flds);
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

