
#pragma once

#include "setup_particles.hxx"

template<>
template<typename FUNC>
void SetupParticles<MparticlesCuda>::setup_particles(Mparticles& mprts, psc* psc,
						     std::vector<uint>& n_prts_by_patch,
						     FUNC func)
{
  auto& mp = mprts.get_as<MparticlesSingle>();
  SetupParticles<MparticlesSingle>::setup_particles(mp, psc, n_prts_by_patch, func);
  mprts.put_as(mp);
}

  
