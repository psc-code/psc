
#pragma once

#include "checks.hxx"

struct ChecksCuda : ChecksBase
{
  void continuity_before_particle_push(psc* psc) { assert(0); }
  void continuity_after_particle_push(psc* psc) { assert(0); }
  void gauss(psc* psc) { assert(0); }

  void continuity_before_particle_push(MparticlesCuda& mprts)
  {}
 
  void continuity_after_particle_push(MparticlesCuda& mprts, MfieldsCuda& mflds)
  {}

  void gauss(MparticlesCuda& mprts, MfieldsCuda& mflds)
  {}
};
