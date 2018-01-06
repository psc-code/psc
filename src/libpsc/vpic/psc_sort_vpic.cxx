
#include "psc_sort_private.h"

#include "psc_particles_vpic.h"
#include "vpic_iface.h"

static void
psc_sort_vpic_run(struct psc_sort *sort, struct psc_mparticles *mprts_base)
{
  struct psc_mparticles *mprts = psc_mparticles_get_as(mprts_base, "vpic", 0);
  struct psc_mparticles_vpic *sub = psc_mparticles_vpic(mprts);
  struct psc *psc = ppsc; // FIXME

  Simulation_sort_mprts(sub->sim, sub->vmprts, psc->timestep);

  psc_mparticles_put_as(mprts, mprts_base, 0);
}

// ----------------------------------------------------------------------
// psc_sort: subclass "vpic"

struct psc_sort_ops psc_sort_vpic_ops = {
  .name                  = "vpic",
  .run                   = psc_sort_vpic_run,
};


