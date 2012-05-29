
#include "psc_push_fields_private.h"

// ----------------------------------------------------------------------
// get_ops

static inline struct psc_push_fields_ops *
get_ops(struct psc_fields *flds)
{
  if (psc_fields_ops(flds) == &psc_fields_single_ops) {
    return &psc_push_fields_single_ops;
  } else if (psc_fields_ops(flds) == &psc_fields_cuda_ops) {
    return &psc_push_fields_cuda_ops;
  } else {
    assert(0);
  }
}

// ----------------------------------------------------------------------

static void
psc_push_fields_mix_push_a_E(struct psc_push_fields *push, struct psc_fields *flds)
{
  struct psc_push_fields_ops *ops = get_ops(flds);
  assert(ops->push_a_E);
  ops->push_a_E(push, flds);
}

static void
psc_push_fields_mix_push_a_H(struct psc_push_fields *push, struct psc_fields *flds)
{
  struct psc_push_fields_ops *ops = get_ops(flds);
  assert(ops->push_a_H);
  ops->push_a_H(push, flds);
}

static void
psc_push_fields_mix_push_b_H(struct psc_push_fields *push, struct psc_fields *flds)
{
  struct psc_push_fields_ops *ops = get_ops(flds);
  assert(ops->push_b_H);
  ops->push_b_H(push, flds);
}

static void
psc_push_fields_mix_push_b_E(struct psc_push_fields *push, struct psc_fields *flds)
{
  struct psc_push_fields_ops *ops = get_ops(flds);
  assert(ops->push_b_E);
  ops->push_b_E(push, flds);
}

// ======================================================================
// psc_push_fields: subclass "mix"

struct psc_push_fields_ops psc_push_fields_mix_ops = {
  .name                  = "mix",
  .push_a_E              = psc_push_fields_mix_push_a_E,
  .push_a_H              = psc_push_fields_mix_push_a_H,
  .push_b_H              = psc_push_fields_mix_push_b_H,
  .push_b_E              = psc_push_fields_mix_push_b_E,
};
