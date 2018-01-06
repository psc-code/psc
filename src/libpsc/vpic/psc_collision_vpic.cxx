
#include "psc_collision_private.h"

#include "psc_particles_vpic.h"
#include "psc_method.h"

#include "vpic_iface.h"

// ----------------------------------------------------------------------

struct psc_collision_vpic {
  Simulation *sim;
};

#define psc_collision_vpic(collision) mrc_to_subobj(collision, struct psc_collision_vpic)

// ----------------------------------------------------------------------
// psc_collision_vpic_setup

static void
psc_collision_vpic_setup(struct psc_collision *collision)
{
  struct psc_collision_vpic *sub = psc_collision_vpic(collision);
  
  psc_method_get_param_ptr(ppsc->method, "sim", (void **) &sub->sim);
}

// ----------------------------------------------------------------------
// psc_collision_vpic_run

static void
psc_collision_vpic_run(struct psc_collision *collision,
		       struct psc_mparticles *mprts_base)
{
  struct psc_collision_vpic *sub = psc_collision_vpic(collision);

  Simulation_collision_run(sub->sim);
}

// ----------------------------------------------------------------------
// psc_collision: subclass "vpic"

struct psc_collision_ops psc_collision_vpic_ops = {
  .name                  = "vpic",
  .size                  = sizeof(struct psc_collision_vpic),
  .setup                 = psc_collision_vpic_setup,
  .run                   = psc_collision_vpic_run,
};

