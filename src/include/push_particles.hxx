
#pragma once

#include "psc_push_particles_private.h"
#include "fields3d.hxx"

// ======================================================================
// PscPushParticlesBase

template<typename S>
struct PscPushParticlesBase
{
  using sub_t = S;

  explicit PscPushParticlesBase(psc_push_particles *pushp)
    : pushp_(pushp),
      sub_(mrc_to_subobj(pushp, sub_t))
  {}

  void operator()(mparticles_base_t mprts, mfields_base_t mflds)
  {
    psc_push_particles_run(pushp_, mprts.mprts(), mflds.mflds());
  }

  psc_push_particles *pushp() { return pushp_; }
  
  sub_t* operator->() { return sub_; }

  sub_t* sub() { return sub_; }

private:
  psc_push_particles *pushp_;
  sub_t *sub_;
};

