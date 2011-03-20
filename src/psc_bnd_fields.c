
#include "psc_bnd_fields_private.h"

// ======================================================================
// forward to subclass

void
psc_bnd_fields_fill_ghosts_H_b(struct psc_bnd_fields *bnd, mfields_base_t *flds)
{
  struct psc_bnd_fields_ops *ops = psc_bnd_fields_ops(bnd);
  assert(ops->fill_ghosts_H_b);
  ops->fill_ghosts_H_b(bnd, flds);
}

// ======================================================================
// psc_bnd_fields_init

static void
psc_bnd_fields_init()
{
  mrc_class_register_subclass(&mrc_class_psc_bnd_fields, &psc_bnd_fields_fortran_ops);
}

// ======================================================================
// psc_bnd_fields class

struct mrc_class_psc_bnd_fields mrc_class_psc_bnd_fields = {
  .name             = "psc_bnd_fields",
  .size             = sizeof(struct psc_bnd_fields),
  .init             = psc_bnd_fields_init,
};

