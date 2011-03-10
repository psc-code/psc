
#include "psc_bnd_private.h"

// ======================================================================
// psc_bnd class

static LIST_HEAD(psc_bnd_subclasses);

void
psc_bnd_register(struct psc_bnd_ops *ops)
{
  list_add_tail(&ops->list, &psc_bnd_subclasses);
}

// ----------------------------------------------------------------------
// psc_bnd_init

static void
psc_bnd_init()
{
  psc_bnd_register_c();
}

struct mrc_class mrc_class_psc_bnd = {
  .name             = "psc_bnd",
  .size             = sizeof(struct psc_bnd),
  .subclasses       = &psc_bnd_subclasses,
  .init             = psc_bnd_init,
};

