
#include "psc_bnd_private.h"

void create_bnd(void);

static void
psc_bnd_c_setup(struct mrc_obj *obj)
{
  create_bnd();
}

// ======================================================================
// psc_bnd: subclass "c"

struct psc_bnd_ops psc_bnd_c_ops = {
  .name                  = "c",
  .size                  = sizeof(struct psc_bnd_c),
  .setup                 = psc_bnd_c_setup,
};
