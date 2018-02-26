
#include "psc_push_particles_private.h"
#include "psc_generic_c.h"

#include "push_particles.hxx"

#include "push.hxx"

struct PushParticlesGenericC : PushParticlesBase
{
  void push_mprts_xyz(struct psc_mparticles *mprts, struct psc_mfields *mflds) override
  { return PscPushParticles_<PushParticles__<Config2ndXYZ>>::push_mprts(nullptr, mprts, mflds); }

  void push_mprts_xy(struct psc_mparticles *mprts, struct psc_mfields *mflds) override
  { return PscPushParticles_<PushParticles__<Config2ndXY>>::push_mprts(nullptr, mprts, mflds); }

  void push_mprts_xz(struct psc_mparticles *mprts, struct psc_mfields *mflds) override
  { return PscPushParticles_<PushParticles__<Config2ndXZ>>::push_mprts(nullptr, mprts, mflds); }

  void push_mprts_yz(struct psc_mparticles *mprts, struct psc_mfields *mflds) override
  { return PscPushParticles_<PushParticles__<Config2ndYZ>>::push_mprts(nullptr, mprts, mflds); }

  void push_mprts_y(struct psc_mparticles *mprts, struct psc_mfields *mflds) override
  { return PscPushParticles_<PushParticles__<Config2ndY>>::push_mprts(nullptr, mprts, mflds); }

  void push_mprts_z(struct psc_mparticles *mprts, struct psc_mfields *mflds) override
  { return PscPushParticles_<PushParticles__<Config2ndZ>>::push_mprts(nullptr, mprts, mflds); }
};

using PushParticlesWrapper_t = PushParticlesWrapper<PushParticlesGenericC>;

// ======================================================================
// psc_push_particles: subclass "generic_c"

struct psc_push_particles_ops_c : psc_push_particles_ops {
  psc_push_particles_ops_c() {
    name                  = "generic_c";
    size                  = PushParticlesWrapper_t::size;
    setup                 = PushParticlesWrapper_t::setup;
    destroy               = PushParticlesWrapper_t::destroy;
    particles_type        = "double";
  }
} psc_push_particles_generic_c_ops;
