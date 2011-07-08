#include "psc_domainwindow_private.h"
#include <mrc_params.h>

#define psc_domainwindow_ops(c) ((struct psc_domainwindow_ops *)((c)->obj.ops))

#define VAR(x) (void *)offsetof(struct psc_domainwindow, x)
static struct param psc_domainwindow_descr[] = {
  { "np"	, VAR(np)	, PARAM_INT3(0,0,0) },
  {}
};
#undef VAR

static void psc_domainwindow_init()
{
  mrc_class_register_subclass(&mrc_class_psc_domainwindow, &psc_movingwindow_z_ops);
}

void psc_domainwindow_timestep(struct psc_domainwindow* this, int timestep, double t)
{
  if (psc_domainwindow_ops(this)->timestep != 0) psc_domainwindow_ops(this)->timestep(this, timestep, t);
}

void _psc_domainwindow_setup(struct psc_domainwindow* this)
{
  bitfield3d_create(&this->activepatches, (unsigned int *) this->np);
  if (psc_domainwindow_ops(this)->setup != 0) psc_domainwindow_ops(this)->setup(this);
}

struct mrc_class_psc_domainwindow mrc_class_psc_domainwindow = {
  .name             = "psc_domainwindow",
  .size             = sizeof(struct psc_domainwindow),
  .param_descr      = psc_domainwindow_descr,
  .setup            = _psc_domainwindow_setup,
  .init             = psc_domainwindow_init,
};
