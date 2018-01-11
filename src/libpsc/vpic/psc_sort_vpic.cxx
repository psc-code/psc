
#include "psc_sort_private.h"

#include "psc_particles_vpic.h"
#include "vpic_iface.h"

static void
psc_sort_vpic_run(struct psc_sort *sort, struct psc_mparticles *mprts_base)
{
  mparticles_vpic_t mprts = mprts_base->get_as<mparticles_vpic_t>();
  struct psc *psc = ppsc; // FIXME

  Simulation_sort_mprts(mprts.sub()->sim, mprts.sub()->vmprts, psc->timestep);

  mprts.put_as(mprts_base);
}

// ----------------------------------------------------------------------
// psc_sort: subclass "vpic"

struct psc_sort_ops_vpic : psc_sort_ops {
  psc_sort_ops_vpic() {
    name                  = "vpic";
    run                   = psc_sort_vpic_run;
  }
} psc_sort_vpic_ops;


