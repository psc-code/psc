
#pragma once

#include "setup_particles.hxx"

#include "psc_particles_cuda.h"
#include "psc_particles_single.h"

template<>
template<typename FUNC>
void SetupParticles<MparticlesCuda<BS144>>::setup_particles(Mparticles& mprts,
							     std::vector<uint>& n_prts_by_patch,
							     FUNC func)
{
  SetupParticles<MparticlesSingle> setup_particles;
  setup_particles.neutralizing_population = neutralizing_population;
  setup_particles.fractional_n_particles_per_cell = fractional_n_particles_per_cell;
  setup_particles.const_num_particles_per_cell = const_num_particles_per_cell;
  setup_particles.initial_momentum_gamma_correction = initial_momentum_gamma_correction;

  auto& mp = mprts.get_as<MparticlesSingle>();
  setup_particles.setup_particles(mp, n_prts_by_patch, func);
  mprts.put_as(mp);
}

// FIXME, missing for BS444

