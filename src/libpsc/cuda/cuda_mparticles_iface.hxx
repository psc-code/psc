
#pragma once

#include "grid.hxx"
#include "particle_simple.hxx"

template<typename BS>
struct cuda_mparticles;

template<typename BS>
struct cuda_mparticles_iface
{
  using cuda_mparticles = cuda_mparticles<BS>;
  using real_t = float;
  
  static cuda_mparticles* new_(const Grid_t& grid);
  static void delete_(cuda_mparticles* cmprts);
  static int get_n_prts(cuda_mparticles* cmprts);
  static std::vector<uint> get_size_all(const cuda_mparticles* cmprts);
  static void inject(cuda_mparticles* cmprts, const std::vector<ParticleSimple<real_t>>& buf,
		     const std::vector<uint>& buf_n_by_patch);
};

