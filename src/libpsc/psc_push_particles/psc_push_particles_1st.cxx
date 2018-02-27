
#include "psc_push_particles_private.h"
#include "psc_push_particles_1st.h"

#include "push_particles.hxx"
#include "push_config.hxx"

#include "psc_particles_as_double.h"
#include "psc_fields_as_c.h"

#include "push_part_common.c"

// ======================================================================
// psc_push_particles: subclass "1st"

template<typename DIM>
using Push = PscPushParticles_<PushParticles__<Config1st<DIM>>>;

struct PushParticles1st : PushParticlesBase
{
  void push_mprts_xz(struct psc_mparticles *mprts, struct psc_mfields *mflds) override
  { return Push<dim_xz>::push_mprts(mprts, mflds); }

  void push_mprts_yz(struct psc_mparticles *mprts, struct psc_mfields *mflds) override
  { return Push<dim_yz>::push_mprts(mprts, mflds); }
};

using PushParticlesWrapper_t = PushParticlesWrapper<PushParticles1st>;

struct psc_push_particles_ops_1st : psc_push_particles_ops {
  psc_push_particles_ops_1st() {
    name                  = "1st";
    size                  = PushParticlesWrapper_t::size;
    setup                 = PushParticlesWrapper_t::setup;
    destroy               = PushParticlesWrapper_t::destroy;
    particles_type        = "double";
  }
} psc_push_particles_1st_ops;

