
#include "push_particles_vpic.hxx"

#include "psc_push_particles_private.h"

// ----------------------------------------------------------------------
// psc_push_particles_vpic_prep

static void
psc_push_particles_vpic_prep(struct psc_push_particles *push,
			     MparticlesBase& mprts_base,
			     MfieldsBase& mflds_base)
{
  PscPushParticles<PushParticlesVpic> pushp(push);
  pushp->prep(mprts_base, mflds_base);
}

// ----------------------------------------------------------------------
// psc_push_particles: subclass "vpic"

using PushParticlesWrapper_t = PushParticlesWrapper<PushParticlesVpic>;

struct psc_push_particles_ops_vpic : psc_push_particles_ops {
  psc_push_particles_ops_vpic() {
    name                  = "vpic";
    size                  = PushParticlesWrapper_t::size;
    setup                 = PushParticlesWrapper_t::setup;
    destroy               = PushParticlesWrapper_t::destroy;
    prep                  = psc_push_particles_vpic_prep;
  }
} psc_push_particles_vpic_ops;

