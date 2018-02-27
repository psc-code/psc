
#ifndef PSC_PUSH_PARTICLES_1VB_H
#define PSC_PUSH_PARTICLES_1VB_H

#include <psc_push_particles_private.h>

#include "push_particles.hxx"
#include "fields.hxx"

#include "../inc_defs.h"
#include "../push_config.hxx"

template<typename C>
struct push_p_ops
{
  using Mparticles = typename C::Mparticles;
  using Mfields = typename C::Mfields;
  using mparticles_t = PscMparticles<Mparticles>;
  using mfields_t = PscMfields<Mfields>;
  
  static void push_mprts(struct psc_mparticles *mprts, struct psc_mfields *mflds_base);
  static void stagger_mprts(struct psc_mparticles *mprts, struct psc_mfields *mflds_base);
  static void push_mprts(Mparticles& mprts, Mfields& mflds);
};

// ======================================================================
// PushParticles_

template<template<class> class PUSH_P_OPS>
class PushParticles_ : public PushParticlesBase
{
  using Self = PushParticles_<PUSH_P_OPS>;
  
public:
  void push_mprts_xyz(struct psc_mparticles *mprts, struct psc_mfields *mflds) override
  { PUSH_P_OPS<dim_xyz>::push_mprts(mprts, mflds); }

#if 0
  void push_mprts_xz(struct psc_mparticles *mprts, struct psc_mfields *mflds) override
  { PUSH_P_OPS<dim_xz>::push_mprts(mprts, mflds_base); }
#endif

  void push_mprts_yz(struct psc_mparticles *mprts, struct psc_mfields *mflds) override
  { PUSH_P_OPS<dim_yz>::push_mprts(mprts, mflds); }

  void push_mprts_1(struct psc_mparticles *mprts, struct psc_mfields *mflds) override
  { PUSH_P_OPS<dim_1>::push_mprts(mprts, mflds); }

  void stagger_mprts_yz(struct psc_mparticles *mprts, struct psc_mfields *mflds) override
  { PUSH_P_OPS<dim_yz>::stagger_mprts(mprts, mflds); }

};


#endif
