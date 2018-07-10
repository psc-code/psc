
#include "bnd_particles_vpic.hxx"

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

