
#include "psc_collision_private.h"

#include "vpic_iface.h"

// ----------------------------------------------------------------------
// psc_collision_vpic_run

static void
psc_collision_vpic_run(struct psc_collision *collision,
		       struct psc_mparticles *mprts_base)
{
  vpic_collision_run();
}

// ----------------------------------------------------------------------
// psc_collision: subclass "vpic"

struct psc_collision_ops psc_collision_vpic_ops = {
  .name                  = "vpic",
  .run                   = psc_collision_vpic_run,
};

