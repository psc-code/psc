
#pragma once

#include "psc_fields_single.h"

#define push_p_ops push_p_ops_1vbec_single
#include "psc_push_particles_1vb.h"

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

struct PushParticles1vbecSingle : PushParticles_<push_p_ops_1vbec_single_>
{
  using Self = PushParticles1vbecSingle;
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

