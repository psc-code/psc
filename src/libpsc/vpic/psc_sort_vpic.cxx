
#include "psc_sort_private.h"

#include "psc_particles_vpic.h"
#include "vpic_iface.h"
#include "sort.hxx"

struct SortVpic : SortBase
{
  using mparticles_t = PscMparticlesVpic;
  
  void operator()(mparticles_t mprts)
  {
    Simulation_sort_mprts(mprts->sim, mprts->vmprts, ppsc->timestep);
  }

  void run(struct psc_mparticles *mprts_base) override
  {
    mparticles_t mprts = mprts_base->get_as<mparticles_t>();
    (*this)(mprts);
    mprts.put_as(mprts_base);
  }
};

// ----------------------------------------------------------------------
// psc_sort: subclass "vpic"

struct psc_sort_ops_vpic : psc_sort_ops {
  using PscSort = PscSortWrapper<SortVpic>;
  psc_sort_ops_vpic() {
    name                  = "vpic";
    size                  = PscSort::size;
    setup                 = PscSort::setup;
    destroy               = PscSort::destroy;
  }
} psc_sort_vpic_ops;


