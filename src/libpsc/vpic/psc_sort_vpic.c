
#include "psc_sort_private.h"

#include "vpic_iface.h"

static void
psc_sort_vpic_run(struct psc_sort *sort, struct psc_mparticles *mprts)
{
  vpic_performance_sort();
}

// ----------------------------------------------------------------------
// psc_sort: subclass "vpic"

struct psc_sort_ops psc_sort_vpic_ops = {
  .name                  = "vpic",
  .run                   = psc_sort_vpic_run,
};


