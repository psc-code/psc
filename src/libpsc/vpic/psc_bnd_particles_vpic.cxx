
#include "bnd_particles.hxx"

#include "vpic_iface.h"

struct BndParticlesVpic : BndParticlesBase
{
  BndParticlesVpic(struct mrc_domain *domain, const Grid_t& grid) {}

  void exchange_particles(MparticlesBase& mprts_base) override {}
};

// ----------------------------------------------------------------------
// psc_bnd_particles: subclass "vpic"

struct psc_bnd_particles_ops_vpic : psc_bnd_particles_ops {
  using Wrapper = PscBndParticlesWrapper<BndParticlesVpic>;
  psc_bnd_particles_ops_vpic() {
    name                    = "vpic";
    size                    = Wrapper::size;
    setup                   = Wrapper::setup;
    destroy                 = Wrapper::destroy;
  }
} psc_bnd_particles_vpic_ops;

