
#ifndef PSC_PUSH_PARTICLES_1VB_H
#define PSC_PUSH_PARTICLES_1VB_H

#include <psc_push_particles_private.h>

#include "push_particles.hxx"
#include "fields.hxx"

#include "../inc_defs.h"
#include "../push_config.hxx"

// ======================================================================
// PushParticles_

template<template<class> class PUSH_P_OPS>
class PushParticles_ : public PushParticlesBase
{
  using Self = PushParticles_<PUSH_P_OPS>;
  
public:
  void push_mprts_xyz(PscMparticlesBase mprts, PscMfieldsBase mflds) override
  { PUSH_P_OPS<dim_xyz>::push_mprts(mprts, mflds); }

#if 0
  void push_mprts_xz(PscMparticlesBase mprts, PscMfieldsBase mflds) override
  { PUSH_P_OPS<dim_xz>::push_mprts(mprts, mflds_base); }
#endif

  void push_mprts_yz(PscMparticlesBase mprts, PscMfieldsBase mflds) override
  { push_yz_.push_mprts(mprts, mflds); }

  void push_mprts_1(PscMparticlesBase mprts, PscMfieldsBase mflds) override
  { PUSH_P_OPS<dim_1>::push_mprts(mprts, mflds); }

  void stagger_mprts_yz(PscMparticlesBase mprts, PscMfieldsBase mflds) override
  { PUSH_P_OPS<dim_yz>::stagger_mprts(mprts, mflds); }

  PUSH_P_OPS<dim_yz> push_yz_;
};


#endif
