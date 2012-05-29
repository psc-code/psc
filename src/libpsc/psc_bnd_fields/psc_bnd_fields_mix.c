
#include "psc_bnd_fields_private.h"

// ----------------------------------------------------------------------
// get_ops

static inline struct psc_bnd_fields_ops *
get_ops(struct psc_fields *flds)
{
  if (psc_fields_ops(flds) == &psc_fields_single_ops) {
    return &psc_bnd_fields_single_ops;
  } else if (psc_fields_ops(flds) == &psc_fields_cuda_ops) {
    return &psc_bnd_fields_cuda_ops;
  } else {
    assert(0);
  }
}

// ----------------------------------------------------------------------

static void
psc_bnd_fields_mix_fill_ghosts_a_E(struct psc_bnd_fields *bnd, struct psc_fields *flds)
{
  struct psc_bnd_fields_ops *ops = get_ops(flds);
  assert(ops->fill_ghosts_a_E);
  ops->fill_ghosts_a_E(bnd, flds);
}

static void
psc_bnd_fields_mix_fill_ghosts_a_H(struct psc_bnd_fields *bnd, struct psc_fields *flds)
{
  struct psc_bnd_fields_ops *ops = get_ops(flds);
  assert(ops->fill_ghosts_a_H);
  ops->fill_ghosts_a_H(bnd, flds);
}

static void
psc_bnd_fields_mix_fill_ghosts_b_E(struct psc_bnd_fields *bnd, struct psc_fields *flds)
{
  struct psc_bnd_fields_ops *ops = get_ops(flds);
  assert(ops->fill_ghosts_b_E);
  ops->fill_ghosts_b_E(bnd, flds);
}

static void
psc_bnd_fields_mix_fill_ghosts_b_H(struct psc_bnd_fields *bnd, struct psc_fields *flds)
{
  struct psc_bnd_fields_ops *ops = get_ops(flds);
  assert(ops->fill_ghosts_b_H);
  ops->fill_ghosts_b_H(bnd, flds);
}

static void
psc_bnd_fields_mix_add_ghosts_J(struct psc_bnd_fields *bnd, struct psc_fields *flds)
{
  struct psc_bnd_fields_ops *ops = get_ops(flds);
  assert(ops->add_ghosts_J);
  ops->add_ghosts_J(bnd, flds);
}

// ======================================================================
// psc_bnd_fields: subclass "mix"

struct psc_bnd_fields_ops psc_bnd_fields_mix_ops = {
  .name                  = "mix",
  .fill_ghosts_a_E       = psc_bnd_fields_mix_fill_ghosts_a_E,
  .fill_ghosts_a_H       = psc_bnd_fields_mix_fill_ghosts_a_H,
  .fill_ghosts_b_E       = psc_bnd_fields_mix_fill_ghosts_b_E,
  .fill_ghosts_b_H       = psc_bnd_fields_mix_fill_ghosts_b_H,
  .add_ghosts_J          = psc_bnd_fields_mix_add_ghosts_J,
};
