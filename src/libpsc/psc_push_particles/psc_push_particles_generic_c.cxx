
#include "psc_push_particles_private.h"

#include "push_particles.hxx"
#include "push_config.hxx"

#include "psc_particles_as_double.h"
#include "psc_fields_as_c.h"

#include "push_part_common.c"

template<typename DIM>
using Push = PscPushParticles_<PushParticles__<Config2nd<DIM>>>;

struct PushParticlesGenericC : PushParticlesBase
{
  void push_mprts_xyz(struct psc_mparticles *mprts, struct psc_mfields *mflds) override
  { return Push<dim_xyz>::push_mprts(mprts, mflds); }

  void push_mprts_xy(struct psc_mparticles *mprts, struct psc_mfields *mflds) override
  { return Push<dim_xy>::push_mprts(mprts, mflds); }

  void push_mprts_xz(struct psc_mparticles *mprts, struct psc_mfields *mflds) override
  { return Push<dim_xz>::push_mprts(mprts, mflds); }

  void push_mprts_yz(struct psc_mparticles *mprts, struct psc_mfields *mflds) override
  { return Push<dim_yz>::push_mprts(mprts, mflds); }

  void push_mprts_y(struct psc_mparticles *mprts, struct psc_mfields *mflds) override
  { return Push<dim_y>::push_mprts(mprts, mflds); }

  void push_mprts_z(struct psc_mparticles *mprts, struct psc_mfields *mflds) override
  { return Push<dim_z>::push_mprts(mprts, mflds); }

  void push_mprts_1(struct psc_mparticles *mprts, struct psc_mfields *mflds) override
  { return Push<dim_1>::push_mprts(mprts, mflds); }
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
