
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
  using Particle = ParticleSimple<real_t>;
  
  static cuda_mparticles* new_(const Grid_t& grid);
  static void delete_(cuda_mparticles* cmprts);
  static int get_n_prts(cuda_mparticles* cmprts);
  static std::vector<uint> get_size_all(const cuda_mparticles* cmprts);
  static void inject(cuda_mparticles* cmprts, const std::vector<Particle>& buf,
		     const std::vector<uint>& buf_n_by_patch);

  static std::vector<uint> get_offsets(const cuda_mparticles* cmprts);
  static std::vector<Particle> get_particles(const cuda_mparticles* cmprts);
  static Particle get_particle(const cuda_mparticles* cmprts, int p, int n);

  static bool check_after_push(cuda_mparticles* cmprts);
  static void dump(const cuda_mparticles* cmprts, const std::string& filename);
};

