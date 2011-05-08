
#include "psc_push_fields_private.h"

#include "psc_bnd.h"
#include "psc_bnd_fields.h"
#include <mrc_profile.h>

static void
_psc_push_fields_create(struct psc_push_fields *push)
{
  push->bnd_fields = psc_bnd_fields_create(psc_push_fields_comm(push));
}

static void
_psc_push_fields_set_from_options(struct psc_push_fields *push)
{
  psc_bnd_fields_set_from_options(push->bnd_fields);
}

static void
_psc_push_fields_setup(struct psc_push_fields *push)
{
  // The cbe version (and possibly cuda) need
  // the subclass (ops) setup to be called at
  // some point. So we do need to have this forwarded
  // to the subclass. --sra
  if(psc_push_fields_ops(push)->setup) {
    psc_push_fields_ops(push)->setup(push);
  }
  psc_bnd_fields_setup(push->bnd_fields);
}

static void
_psc_push_fields_view(struct psc_push_fields *push)
{
  psc_bnd_fields_view(push->bnd_fields);
}

static void
_psc_push_fields_destroy(struct psc_push_fields *push)
{
  psc_bnd_fields_destroy(push->bnd_fields);
  if(psc_push_fields_ops(push)->destroy) {
    psc_push_fields_ops(push)->destroy(push);
  }
}

struct psc_bnd_fields *
psc_push_fields_get_bnd_fields(struct psc_push_fields *push)
{
  return push->bnd_fields;
}

// ======================================================================
// forward to subclass

static inline void
psc_push_fields_push_a_E(struct psc_push_fields *push, mfields_base_t *flds)
{
  struct psc_push_fields_ops *ops = psc_push_fields_ops(push);
  static int pr;
  if (!pr) {
    pr = prof_register("push_fields_a_E", 1., 0, 0);
  }
  assert(ops->push_a_E);
  
  prof_start(pr);
  ops->push_a_E(push, flds);
  prof_stop(pr);
}

static inline void
psc_push_fields_push_a_H(struct psc_push_fields *push, mfields_base_t *flds)
{
  struct psc_push_fields_ops *ops = psc_push_fields_ops(push);
  static int pr;
  if (!pr) {
    pr = prof_register("push_fields_a_H", 1., 0, 0);
  }
  assert(ops->push_a_H);
  
  prof_start(pr);
  ops->push_a_H(push, flds);
  prof_stop(pr);
}

static inline void
psc_push_fields_push_b_E(struct psc_push_fields *push, mfields_base_t *flds)
{
  struct psc_push_fields_ops *ops = psc_push_fields_ops(push);
  static int pr;
  if (!pr) {
    pr = prof_register("push_fields_b_E", 1., 0, 0);
  }
  assert(ops->push_b_E);
  
  prof_start(pr);
  ops->push_b_E(push, flds);
  prof_stop(pr);
}

static inline void
psc_push_fields_push_b_H(struct psc_push_fields *push, mfields_base_t *flds)
{
  struct psc_push_fields_ops *ops = psc_push_fields_ops(push);
  static int pr;
  if (!pr) {
    pr = prof_register("push_fields_b_H", 1., 0, 0);
  }
  assert(ops->push_b_H);
  
  prof_start(pr);
  ops->push_b_H(push, flds);
  prof_stop(pr);
}

void
psc_push_fields_step_a(struct psc_push_fields *push, mfields_base_t *flds)
{
  if (ppsc->domain.use_pml) {
    // FIXME, pml routines sehould be split into E, H push + ghost points, too
    // pml could become a separate push_fields subclass
    struct psc_push_fields_ops *ops = psc_push_fields_ops(push);
    assert(ops->pml_a);
    ops->pml_a(push, flds);
  } else {
    psc_push_fields_push_a_E(push, flds);
    psc_bnd_fill_ghosts(ppsc->bnd, flds, EX, EX + 3);

    psc_push_fields_push_a_H(push, flds);
    psc_bnd_fill_ghosts(ppsc->bnd, flds, HX, HX + 3);
  }
}

void
psc_push_fields_step_b(struct psc_push_fields *push, mfields_base_t *flds)
{
  if (ppsc->domain.use_pml) {
    struct psc_push_fields_ops *ops = psc_push_fields_ops(push);
    assert(ops->pml_b);
    ops->pml_b(push, flds);
  } else {
    psc_push_fields_push_b_H(push, flds);
    psc_bnd_fields_fill_ghosts_b_H(push->bnd_fields, flds);
    psc_bnd_fill_ghosts(ppsc->bnd, flds, HX, HX + 3);

    psc_push_fields_push_b_E(push, flds);
    psc_bnd_fill_ghosts(ppsc->bnd, flds, EX, EX + 3);
  }
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
  .create           = _psc_push_fields_create,
  .set_from_options = _psc_push_fields_set_from_options,
  .setup            = _psc_push_fields_setup,
  .view             = _psc_push_fields_view,
  .destroy          = _psc_push_fields_destroy,
};

