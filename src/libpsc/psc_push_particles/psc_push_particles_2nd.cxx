
#include "psc_push_particles_private.h"

#include "push_particles.hxx"
#include "push_config.hxx"

#include "push_part_common.c"

// ======================================================================
// psc_push_particles: subclass "2nd_double"
//
// 2nd order, Esirkepov, like the old generic_c/fortran, but on "double"
// particles

struct PushParticles2ndDouble : PushParticlesBase
{
  void push_mprts_yz(struct psc_mparticles *mprts, struct psc_mfields *mflds) override
  { return PscPushParticles_<PushParticles__<Config2ndDoubleYZ>>::push_mprts(mprts, mflds); }
};

using PushParticlesWrapper_t = PushParticlesWrapper<PushParticles2ndDouble>;

struct psc_push_particles_ops_2nd_double : psc_push_particles_ops {
  psc_push_particles_ops_2nd_double() {
    name                  = "2nd_double";
    size                  = PushParticlesWrapper_t::size;
    setup                 = PushParticlesWrapper_t::setup;
    destroy               = PushParticlesWrapper_t::destroy;
  }
} psc_push_particles_2nd_double_ops;

