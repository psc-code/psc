
#include "psc_sort_private.h"

#include "psc_particles_vpic.h"
#include "vpic_iface.h"
#include "sort.hxx"

struct SortVpic : SortCRTP<SortVpic, MparticlesVpic>
{
  using Base = SortCRTP<SortVpic, MparticlesVpic>;
  using Mparticles = MparticlesVpic;

  using Base::Base;
  
  void operator()(Mparticles& mprts)
  {
    Simulation_sort_mprts(mprts.sim, mprts.vmprts, ppsc->timestep);
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


