
#include "psc_bnd_private.h"

// ======================================================================
// psc_bnd_init

static void
psc_bnd_init()
{
  psc_bnd_register_c();
}

// ======================================================================
// psc_bnd class

struct mrc_class mrc_class_psc_bnd = {
  .name             = "psc_bnd",
  .size             = sizeof(struct psc_bnd),
  .init             = psc_bnd_init,
};

