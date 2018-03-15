
#pragma once

#include "psc_fields_single.h"
#include "psc_particles_single.h"

#include "psc_push_particles_1vb.h"

#include "../1vb.c"

template<typename dim>
using push_p_ops_1vbec_double = PscPushParticles_<PushParticles1vb<Config1vbecDouble<dim>>>;

using PushParticles_t = PushParticles_<push_p_ops_1vbec_double>;


template<typename dim>
using push_p_ops_1vbec_single_ = PscPushParticles_<PushParticles1vb<Config1vbecSingle<dim>>>;

// FIXME, special hack... for xyz_xz

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

