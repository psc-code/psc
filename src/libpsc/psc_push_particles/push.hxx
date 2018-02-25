
#pragma once

#include "push_config.hxx"

template<typename C>
struct PushParticles__;

template<typename PushParticles_t>
struct PscPushParticles_
{
  using Mparticles = MparticlesDouble; // FIXME, horrible...
  using Mfields = MfieldsC;
  
  static void push_mprts(struct psc_push_particles *push,
			 struct psc_mparticles *mprts,
			 struct psc_mfields *mflds_base);
  static void push_mprts(Mparticles& mprts, Mfields& mflds);
};

