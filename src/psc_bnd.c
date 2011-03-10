
#include "psc_bnd_private.h"

// ======================================================================
// psc_bnd_init

static void
psc_bnd_init()
{
  mrc_class_register_subclass(&mrc_class_psc_bnd, &psc_bnd_c_ops);
}

// ======================================================================
// psc_bnd class

struct mrc_class mrc_class_psc_bnd = {
  .name             = "psc_bnd",
  .size             = sizeof(struct psc_bnd),
  .init             = psc_bnd_init,
};

