
#pragma once

#include "push_config.hxx"

template<typename C>
struct PushParticles_;

template<typename PushParticles_t>
struct PscPushParticles_
{
  static void push_mprts(struct psc_push_particles *push,
			 struct psc_mparticles *mprts,
			 struct psc_mfields *mflds_base);
};

