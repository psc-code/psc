
#pragma once

#include "psc_fields_single.h"
#include "psc_particles_single.h"

#include "psc_push_particles_1vb.h"

#include "../1vb.c"

template<typename Mparticles, typename Mfields, typename dim>
using push_p_ops_1vbec = PscPushParticles_<PushParticles1vb<Config1vbec<Mparticles, Mfields, dim>>>;

template<typename dim>
using push_p_ops_1vbec_double = push_p_ops_1vbec<MparticlesDouble, MfieldsC, dim>;

using PushParticles1vbecDouble_t = PushParticles_<push_p_ops_1vbec_double>;

template<typename dim>
using push_p_ops_1vbec_single = push_p_ops_1vbec<MparticlesSingle, MfieldsSingle, dim>;

// FIXME, special hack... for xyz_xz

using PushParticles1vbecSingle_t = PushParticles_<push_p_ops_1vbec_single>;

#if 0
struct PushParticles1vbecSingle : PushParticles_<push_p_ops_1vbec_single_>
{
  template<typename dim_t>
  using ops_t = push_p_ops_1vbec_single_<dim_t>;

  using Self = PushParticles1vbecSingle;
  using Base = PushParticles_<push_p_ops_1vbec_single_>;

  using Base::push_mprts_xz;

  void push_mprts_xz(PscMparticlesBase mprts, PscMfieldsBase mflds) override
  {
    PscPushParticles_<PushParticles1vb<Config1vbecSingleXZ>>::push_mprts(mprts, mflds);
  }
};
#endif

