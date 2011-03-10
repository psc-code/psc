
#include "psc_bnd_private.h"

struct psc_bnd_c {
};

// ======================================================================
// psc_bnd: subclass "c"

static struct psc_bnd_ops psc_bnd_c_ops = {
  .name                  = "c",
  .size                  = sizeof(struct psc_bnd_c),
};

void
psc_bnd_register_c()
{
  mrc_class_register_subclass(&mrc_class_psc_bnd, &psc_bnd_c_ops);
}
