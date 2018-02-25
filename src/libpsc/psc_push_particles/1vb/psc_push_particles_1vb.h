
#ifndef PSC_PUSH_PARTICLES_1VB_H
#define PSC_PUSH_PARTICLES_1VB_H

#include <psc_push_particles_private.h>

#include "push_particles.hxx"
#include "fields.hxx"

#include "../inc_defs.h"

template<typename MP, typename MF, typename D, typename O, typename J,
	 typename OPT_EXT = opt_ext_none>
struct push_p_config
{
  using mparticles_t = MP;
  using mfields_t = MF;
  using dim_t = D;
  using order_t = O;
  using calcj_t = J;
  using ext_t = OPT_EXT;
};

template<typename C>
struct push_p_ops
{
  using mparticles_t = typename C::mparticles_t;
  using mfields_t = typename C::mfields_t;
  using Mparticles = typename mparticles_t::sub_t;
  using Mfields = typename mfields_t::sub_t;
  
  static void push_mprts(struct psc_push_particles *push, struct psc_mparticles *mprts,
			 struct psc_mfields *mflds_base);
  static void stagger_mprts(struct psc_push_particles *push, struct psc_mparticles *mprts,
			    struct psc_mfields *mflds_base);
  static void push_mprts(Mparticles& mprts, Mfields& mflds);
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
