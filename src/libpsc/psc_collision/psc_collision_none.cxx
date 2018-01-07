
#include "psc_collision_private.h"

#include <mrc_profile.h>

// ----------------------------------------------------------------------
// psc_collision_none_run

static void
psc_collision_none_run(struct psc_collision *collision,
		       struct psc_mparticles *mprts_base)
{
}

// ======================================================================
// psc_collision: subclass "none"

struct psc_collision_ops psc_collision_none_ops = {
  .name                  = "none",
  .run                   = psc_collision_none_run,
};
