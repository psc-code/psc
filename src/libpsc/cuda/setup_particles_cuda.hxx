
#pragma once

#include "setup_particles.hxx"

template<>
template<typename FUNC>
void SetupParticles<MparticlesCuda<BS144>>::setup_particles(Mparticles& mprts, psc* psc,
							    std::vector<uint>& n_prts_by_patch,
							    FUNC func)
{
  auto& mp = mprts.get_as<MparticlesSingle>();
  SetupParticles<MparticlesSingle>::setup_particles(mp, psc, n_prts_by_patch, func);
  mprts.put_as(mp);
}

template<>
template<typename FUNC>
void SetupParticles<MparticlesCuda<BS144>>::setup_particles(MparticlesCuda<BS144>& mprts,
							    std::vector<uint>& n_prts_by_patch,
							    FUNC func)
{
  using particle_t = typename Mparticles::particle_t;

  std::vector<cuda_mparticles_prt> buf;
  for (int p = 0; p < mprts.n_patches(); p++) {
    for (int n = 0; n < n_prts_by_patch[p]; n++) {
      particle_t prt = func(p, n);
      cuda_mparticles_prt cprt;
      cprt.xi[0] = prt.xi;
      cprt.xi[1] = prt.yi;
      cprt.xi[2] = prt.zi;
      cprt.pxi[0] = prt.pxi;
      cprt.pxi[1] = prt.pyi;
      cprt.pxi[2] = prt.pzi;
      cprt.kind = prt.kind_;
      cprt.qni_wni = prt.qni_wni_;
      buf.push_back(cprt);
    }
  }
    
  mprts.reserve_all(n_prts_by_patch.data());
  mprts.inject_buf(buf.data(), n_prts_by_patch.data());
}
