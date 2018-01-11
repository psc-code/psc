
#include "psc_bnd_fields_private.h"
#include "psc_fields_as_c.h"

// ======================================================================
// forward to subclass

void
psc_bnd_fields_fill_ghosts_E(struct psc_bnd_fields *bnd, struct psc_mfields *mflds)
{
  struct psc_bnd_fields_ops *ops = psc_bnd_fields_ops(bnd);
  if (ops->fill_ghosts_E) {
    ops->fill_ghosts_E(bnd, mflds);
  }
}

void
psc_bnd_fields_fill_ghosts_H(struct psc_bnd_fields *bnd, struct psc_mfields *mflds)
{
  struct psc_bnd_fields_ops *ops = psc_bnd_fields_ops(bnd);
  if (ops->fill_ghosts_H) {
    ops->fill_ghosts_H(bnd, mflds);
  }
}

void
psc_bnd_fields_add_ghosts_J(struct psc_bnd_fields *bnd, struct psc_mfields *mflds)
{
  struct psc_bnd_fields_ops *ops = psc_bnd_fields_ops(bnd);
  if (ops->add_ghosts_J) {
    ops->add_ghosts_J(bnd, mflds);
  }
}

// ======================================================================
// psc_bnd_fields_init

extern struct psc_bnd_fields_ops psc_bnd_fields_fortran_ops;
extern struct psc_bnd_fields_ops psc_bnd_fields_none_ops;
extern struct psc_bnd_fields_ops psc_bnd_fields_c_ops;
extern struct psc_bnd_fields_ops psc_bnd_fields_single_ops;
extern struct psc_bnd_fields_ops psc_bnd_fields_cuda_ops;
extern struct psc_bnd_fields_ops psc_bnd_fields_vpic_ops;

static void
psc_bnd_fields_init()
{
  mrc_class_register_subclass(&mrc_class_psc_bnd_fields, &psc_bnd_fields_c_ops);
  mrc_class_register_subclass(&mrc_class_psc_bnd_fields, &psc_bnd_fields_single_ops);
  mrc_class_register_subclass(&mrc_class_psc_bnd_fields, &psc_bnd_fields_none_ops);
#ifdef USE_FORTRAN
  mrc_class_register_subclass(&mrc_class_psc_bnd_fields, &psc_bnd_fields_fortran_ops);
#endif
#ifdef USE_CUDA
  mrc_class_register_subclass(&mrc_class_psc_bnd_fields, &psc_bnd_fields_cuda_ops);
#endif
#ifdef USE_VPIC
  mrc_class_register_subclass(&mrc_class_psc_bnd_fields, &psc_bnd_fields_vpic_ops);
#endif
}

// ======================================================================
// psc_bnd_fields class

struct mrc_class_psc_bnd_fields_ : mrc_class_psc_bnd_fields {
  mrc_class_psc_bnd_fields_() {
    name             = "psc_bnd_fields";
    size             = sizeof(struct psc_bnd_fields);
    init             = psc_bnd_fields_init;
  }
} mrc_class_psc_bnd_fields;

