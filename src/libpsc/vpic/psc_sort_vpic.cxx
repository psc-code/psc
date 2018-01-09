
#include "psc_sort_private.h"

#include "psc_particles_vpic.h"
#include "vpic_iface.h"

static void
psc_sort_vpic_run(struct psc_sort *sort, struct psc_mparticles *mprts_base)
{
  mparticles_vpic_t mprts = mprts_base->get_as<mparticles_vpic_t>();
  struct psc_mparticles_vpic *sub = psc_mparticles_vpic(mprts.mprts());
  struct psc *psc = ppsc; // FIXME

  Simulation_sort_mprts(sub->sim, sub->vmprts, psc->timestep);

  mprts.put_as(mprts_base);
}

// ----------------------------------------------------------------------
// psc_sort: subclass "vpic"

struct psc_sort_ops psc_sort_vpic_ops = {
  .name                  = "vpic",
  .run                   = psc_sort_vpic_run,
};


