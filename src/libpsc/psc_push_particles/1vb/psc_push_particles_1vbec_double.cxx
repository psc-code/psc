
#include "psc_push_particles_private.h"

#include "psc_particles_as_double.h"
#include "psc_fields_as_c.h"
#include "psc_push_particles_1vb.h"
#include "push_particles.hxx"

// ======================================================================
// psc_push_particles: subclass "1vbec_double"

template<typename dim_t>
using push_p_ops_1vbec_double = push_p_ops<push_p_config<mfields_c_t, dim_t>>;

// ======================================================================
// PushParticles1vbecDouble

template<template<class> class PUSH_P_OPS>
class PushParticles_ : public PushParticlesBase
{
  using Self = PushParticles_<PUSH_P_OPS>;
  
public:
  void push_mprts_xyz(struct psc_mparticles *mprts, struct psc_mfields *mflds_base) override
  { PUSH_P_OPS<dim_xyz>::push_mprts(nullptr, mprts, mflds_base); }

  void push_mprts_yz(struct psc_mparticles *mprts, struct psc_mfields *mflds_base) override
  { PUSH_P_OPS<dim_yz>::push_mprts(nullptr, mprts, mflds_base); }

  void push_mprts_1(struct psc_mparticles *mprts, struct psc_mfields *mflds_base) override
  { PUSH_P_OPS<dim_1>::push_mprts(nullptr, mprts, mflds_base); }

  static void setup(struct psc_push_particles *push)
  {
    PscPushParticles<Self> pushp(push);
    new(pushp.sub()) Self;
  }

  static void destroy(struct psc_push_particles *push)
  {
    PscPushParticles<Self> pushp(push);
    pushp.sub()->~Self();
  }

  static void push_mprts_xyz(struct psc_push_particles *push,
			     struct psc_mparticles *mprts, struct psc_mfields *mflds_base)
  {
    PscPushParticles<Self> pushp(push);
    pushp->push_mprts_xyz(mprts, mflds_base);
  }
  
  static void push_mprts_yz(struct psc_push_particles *push,
			    struct psc_mparticles *mprts, struct psc_mfields *mflds_base)
  {
    PscPushParticles<Self> pushp(push);
    pushp->push_mprts_yz(mprts, mflds_base);
  }
  
  static void push_mprts_1(struct psc_push_particles *push,
			   struct psc_mparticles *mprts, struct psc_mfields *mflds_base)
  {
    PscPushParticles<Self> pushp(push);
    pushp->push_mprts_1(mprts, mflds_base);
  }

};

using PushParticles_t = PushParticles_<push_p_ops_1vbec_double>;
  
struct psc_push_particles_ops_1vbec_double : psc_push_particles_ops {
  psc_push_particles_ops_1vbec_double() {
    name                  = "1vbec_double";
    size                  = sizeof(PushParticles_t);
    setup                 = PushParticles_t::setup;
    destroy               = PushParticles_t::destroy;
    push_mprts_xyz        = PushParticles_t::push_mprts_xyz;
    push_mprts_yz         = PushParticles_t::push_mprts_yz;
    push_mprts_1          = PushParticles_t::push_mprts_1;
    //stagger_mprts_1      = push_p_ops_1vbec_double<dim_1>::stagger_mprts;
    particles_type        = PARTICLE_TYPE;
  }
} psc_push_particles_1vbec_double_ops;

