
#pragma once

#include "../libpsc/psc_push_particles/push_config.hxx"

#include "push_part_common.c"

template<typename dim>
using PushGenericC = PscPushParticles_<PushParticles__<Config2ndDouble<dim>>>;

struct PushParticlesGenericC : PushParticlesBase
{
  void push_mprts_xyz(MparticlesBase& mprts, MfieldsBase& mflds) override
  { return PushGenericC<dim_xyz>::push_mprts(mprts, mflds); }

  void push_mprts_xy(MparticlesBase& mprts, MfieldsBase& mflds) override
  { return PushGenericC<dim_xy>::push_mprts(mprts, mflds); }

  void push_mprts_xz(MparticlesBase& mprts, MfieldsBase& mflds) override
  { return PushGenericC<dim_xz>::push_mprts(mprts, mflds); }

  void push_mprts_yz(MparticlesBase& mprts, MfieldsBase& mflds) override
  { return PushGenericC<dim_yz>::push_mprts(mprts, mflds); }

  void push_mprts_y(MparticlesBase& mprts, MfieldsBase& mflds) override
  { return PushGenericC<dim_y>::push_mprts(mprts, mflds); }

  void push_mprts_z(MparticlesBase& mprts, MfieldsBase& mflds) override
  { return PushGenericC<dim_z>::push_mprts(mprts, mflds); }

  void push_mprts_1(MparticlesBase& mprts, MfieldsBase& mflds) override
  { return PushGenericC<dim_1>::push_mprts(mprts, mflds); }

  PushGenericC<dim_yz> push_yz_;
};


