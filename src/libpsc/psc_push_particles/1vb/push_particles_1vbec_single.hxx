
#pragma once

#include "psc_fields_single.h"
#include "psc_particles_single.h"

#include "psc_push_particles_1vb.h"

template<typename dim>
using push_p_ops_1vbec_single_ = push_p_ops<Config1vbecDouble<dim>>;

// FIXME, special hack... for xyz_xz

struct PushParticles1vbecSingle : PushParticles_<push_p_ops_1vbec_single_>
{
  template<typename dim_t>
  using ops_t = push_p_ops_1vbec_single_<dim_t>;

  using Self = PushParticles1vbecSingle;
  using Base = PushParticles_<push_p_ops_1vbec_single_>;

  using Base::push_mprts_xz;

  void push_mprts_xz(struct psc_mparticles *mprts, struct psc_mfields *mflds_base) override
  {
    push_p_ops<Config1vbecSingleXZ>::push_mprts(nullptr, mprts, mflds_base);
  }

  static void push_mprts_xz(struct psc_push_particles *push,
			    struct psc_mparticles *mprts, struct psc_mfields *mflds_base)
  {
    PscPushParticles<Self> pushp(push);
    pushp->push_mprts_xz(mprts, mflds_base);
  }

};

