
#include "psc_push_particles_private.h"
#include "psc_push_particles_1st.h"

#include "push_particles.hxx"

// ======================================================================
// psc_push_particles: subclass "1st"

struct PushParticles1st : PushParticlesBase
{
  void push_mprts_xz(struct psc_mparticles *mprts, struct psc_mfields *mflds) override
  { return psc_push_particles_push_mprts_1st_xz(nullptr, mprts, mflds); }

  void push_mprts_yz(struct psc_mparticles *mprts, struct psc_mfields *mflds) override
  { return psc_push_particles_push_mprts_1st_yz(nullptr, mprts, mflds); }
};

using PushParticlesWrapper_t = PushParticlesWrapper<PushParticles1st>;

struct psc_push_particles_ops_1st : psc_push_particles_ops {
  psc_push_particles_ops_1st() {
    name                  = "1st";
    size                  = PushParticlesWrapper_t::size;
    setup                 = PushParticlesWrapper_t::setup;
    destroy               = PushParticlesWrapper_t::destroy;
    push_mprts_xz         = PushParticlesWrapper_t::push_mprts_xz;
    push_mprts_yz         = PushParticlesWrapper_t::push_mprts_yz;
    particles_type        = "double";
  }
} psc_push_particles_1st_ops;

