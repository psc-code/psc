
#include "psc_sort_private.h"

#include <mrc_profile.h>

// ----------------------------------------------------------------------
// psc_sort_none_run

static void
psc_sort_none_run(struct psc_sort *sort,
		  mparticles_base_t *particles_base)
{
}

// ======================================================================
// psc_sort: subclass "none"

struct psc_sort_ops psc_sort_none_ops = {
  .name                  = "none",
  .run                   = psc_sort_none_run,
};
