
#pragma once

#include "push_part_common.c"

template<typename dim>
using PushGenericC = PscPushParticles_<PushParticles__<Config2nd<dim>>>;

struct PushParticlesGenericC : PushParticlesBase
{
  void push_mprts_xyz(PscMparticlesBase mprts, PscMfieldsBase mflds) override
  { return PushGenericC<dim_xyz>::push_mprts(mprts, mflds); }

  void push_mprts_xy(PscMparticlesBase mprts, PscMfieldsBase mflds) override
  { return PushGenericC<dim_xy>::push_mprts(mprts, mflds); }

  void push_mprts_xz(PscMparticlesBase mprts, PscMfieldsBase mflds) override
  { return PushGenericC<dim_xz>::push_mprts(mprts, mflds); }

  void push_mprts_yz(PscMparticlesBase mprts, PscMfieldsBase mflds) override
  { return PushGenericC<dim_yz>::push_mprts(mprts, mflds); }

  void push_mprts_y(PscMparticlesBase mprts, PscMfieldsBase mflds) override
  { return PushGenericC<dim_y>::push_mprts(mprts, mflds); }

  void push_mprts_z(PscMparticlesBase mprts, PscMfieldsBase mflds) override
  { return PushGenericC<dim_z>::push_mprts(mprts, mflds); }

  void push_mprts_1(PscMparticlesBase mprts, PscMfieldsBase mflds) override
  { return PushGenericC<dim_1>::push_mprts(mprts, mflds); }

  PushGenericC<dim_yz> push_yz_;
};


