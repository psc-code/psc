
#ifndef PSC_PUSH_PARTICLES_1VB_H
#define PSC_PUSH_PARTICLES_1VB_H

#include <psc_push_particles_private.h>

#include "push_particles.hxx"
#include "fields.hxx"

#include "../inc_defs.h"

template<typename MF, typename D, typename O, typename J>
struct push_p_config
{
  using mfields_t = MF;
  using dim_t = D;
  using order_t = O;
  using calcj_t = J;
};

template<typename C>
struct push_p_ops
{
  static void push_mprts(struct psc_push_particles *push, struct psc_mparticles *mprts,
			 struct psc_mfields *mflds_base);
  static void stagger_mprts(struct psc_push_particles *push, struct psc_mparticles *mprts,
			    struct psc_mfields *mflds_base);
};

// ======================================================================
// PushParticles_

template<template<class> class PUSH_P_OPS>
class PushParticles_ : public PushParticlesBase
{
  using Self = PushParticles_<PUSH_P_OPS>;
  
public:
  void push_mprts_xyz(struct psc_mparticles *mprts, struct psc_mfields *mflds_base) override
  { PUSH_P_OPS<dim_xyz>::push_mprts(nullptr, mprts, mflds_base); }

#if 0
  void push_mprts_xz(struct psc_mparticles *mprts, struct psc_mfields *mflds_base) override
  { PUSH_P_OPS<dim_xz>::push_mprts(nullptr, mprts, mflds_base); }
#endif

  void push_mprts_yz(struct psc_mparticles *mprts, struct psc_mfields *mflds_base) override
  { PUSH_P_OPS<dim_yz>::push_mprts(nullptr, mprts, mflds_base); }

  void push_mprts_1(struct psc_mparticles *mprts, struct psc_mfields *mflds_base) override
  { PUSH_P_OPS<dim_1>::push_mprts(nullptr, mprts, mflds_base); }

  void stagger_mprts_yz(struct psc_mparticles *mprts, struct psc_mfields *mflds_base) override
  { PUSH_P_OPS<dim_yz>::stagger_mprts(nullptr, mprts, mflds_base); }

};


#endif
