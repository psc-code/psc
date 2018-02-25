
#pragma once

#include "psc_fields_single.h"
#include "psc_particles_single.h"

#include "psc_push_particles_1vb.h"

template<typename dim_t>
using push_p_ops_1vbec_single_ = push_p_ops<push_p_config<mparticles_single_t, PscMfieldsSingle, dim_t, opt_order_1st, opt_calcj_1vb_var1>>;

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
    push_p_ops<push_p_config<mparticles_single_t, PscMfieldsSingle, dim_xyz, opt_order_1st, opt_calcj_1vb_split>>::push_mprts(nullptr, mprts, mflds_base);
  }

  static void push_mprts_xz(struct psc_push_particles *push,
			    struct psc_mparticles *mprts, struct psc_mfields *mflds_base)
  {
    PscPushParticles<Self> pushp(push);
    pushp->push_mprts_xz(mprts, mflds_base);
  }

};

