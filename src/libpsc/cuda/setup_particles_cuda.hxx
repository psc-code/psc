
#pragma once

#include "setup_particles.hxx"

#include "psc_particles_cuda.h"
#include "psc_particles_single.h"

template<>
template<typename FUNC>
void SetupParticles<MparticlesCuda<BS144>>::setup_particles2(Mparticles& mprts,
							     std::vector<uint>& n_prts_by_patch,
							     FUNC func)
{
  auto& mp = mprts.get_as<MparticlesSingle>();
  assert(0);
  //SetupParticles<MparticlesSingle>::setup_particles(mp, psc, n_prts_by_patch, func);
  mprts.put_as(mp);
}

template<>
template<typename FUNC>
void SetupParticles<MparticlesCuda<BS144>>::setup_particles(MparticlesCuda<BS144>& mprts,
							    std::vector<uint>& n_prts_by_patch,
							    FUNC func)
{
  using particle_t = typename Mparticles::particle_t;

  auto& grid = mprts.grid();

  std::vector<cuda_mparticles_prt> buf;
  for (int p = 0; p < mprts.n_patches(); p++) {
    for (int n = 0; n < n_prts_by_patch[p]; n++) {
      particle_t prt = func(p, n);
      cuda_mparticles_prt cprt{prt.x, prt.p, prt.w, prt.kind};
      buf.push_back(cprt);
    }
  }
    
  mprts.reserve_all(n_prts_by_patch.data());
  mprts.inject_buf(buf.data(), n_prts_by_patch.data());
}

// FIXME, exactly duplicated for BS444

template<>
template<typename FUNC>
void SetupParticles<MparticlesCuda<BS444>>::setup_particles2(Mparticles& mprts,
							    std::vector<uint>& n_prts_by_patch,
							    FUNC func)
{
  auto& mp = mprts.get_as<MparticlesSingle>();
  assert(0);
  //SetupParticles<MparticlesSingle>::setup_particles(mp, psc, n_prts_by_patch, func);
  mprts.put_as(mp);
}

template<>
template<typename FUNC>
void SetupParticles<MparticlesCuda<BS444>>::setup_particles(MparticlesCuda<BS444>& mprts,
							    std::vector<uint>& n_prts_by_patch,
							    FUNC func)
{
  using particle_t = typename Mparticles::particle_t;

  auto& grid = mprts.grid();

  std::vector<cuda_mparticles_prt> buf;
  for (int p = 0; p < mprts.n_patches(); p++) {
    for (int n = 0; n < n_prts_by_patch[p]; n++) {
      particle_t prt = func(p, n);
      cuda_mparticles_prt cprt{prt.x, prt.p, prt.w, prt.kind};
      buf.push_back(cprt);
    }
  }
    
  mprts.reserve_all(n_prts_by_patch.data());
  mprts.inject_buf(buf.data(), n_prts_by_patch.data());
}
