
#include "psc_push_particles_private.h"

#include "psc_particles_as_single.h"
#include "psc_fields_single.h"

#define push_p_ops push_p_ops_1vbec_single
#include "psc_push_particles_1vb.h"

// ======================================================================
// psc_push_particles: subclass "1vbec_single"

template<typename dim_t>
using push_p_ops_1vbec_single_ = push_p_ops_1vbec_single<push_p_config<mfields_single_t, dim_t>>;

// FIXME, special hack... for xyz_xz
template<typename C>
struct push_p_ops_1vbec_single_xz
{
  static void push_mprts(struct psc_push_particles *push, struct psc_mparticles *mprts,
			 struct psc_mfields *mflds_base);
  static void stagger_mprts(struct psc_push_particles *push, struct psc_mparticles *mprts,
			    struct psc_mfields *mflds_base);
};

struct PushParticles_t : PushParticles_<push_p_ops_1vbec_single_>
{
  using Self = PushParticles_t;
  using Base = PushParticles_<push_p_ops_1vbec_single_>;

  using Base::push_mprts_xz;

  void push_mprts_xz(struct psc_mparticles *mprts, struct psc_mfields *mflds_base) override
  {
    push_p_ops_1vbec_single_xz<push_p_config<mfields_single_t, dim_xyz>>::push_mprts(nullptr, mprts, mflds_base);
  }

  static void push_mprts_xz(struct psc_push_particles *push,
			    struct psc_mparticles *mprts, struct psc_mfields *mflds_base)
  {
    PscPushParticles<Self> pushp(push);
    pushp->push_mprts_xz(mprts, mflds_base);
  }
};

using PushParticlesWrapper_t = PushParticlesWrapper<PushParticles_t>;

struct psc_push_particles_ops_1vbec_single : psc_push_particles_ops {
  psc_push_particles_ops_1vbec_single() {
    name                  = "1vbec_single";
    size                  = PushParticlesWrapper_t::size;
    setup                 = PushParticlesWrapper_t::setup;
    destroy               = PushParticlesWrapper_t::destroy;
    particles_type        = PARTICLE_TYPE;
  }
} psc_push_particles_1vbec_single_ops;

