
#pragma once

#include "grid.hxx"
#include "particle_simple.hxx"

template<typename BS>
struct cuda_mparticles;

// ======================================================================
// cuda_mparticles_iface
//
// This is quite a bit of boilerplate code, with the only purpose being to have
// a separation between nvcc-compiled code and CXX compiled code

template<typename BS>
struct cuda_mparticles_iface
{
  using CudaMparticles = cuda_mparticles<BS>;
  using real_t = float;
  using Particle = ParticleSimple<real_t>;
  
  static CudaMparticles* new_(const Grid_t& grid);
  static void delete_(CudaMparticles* cmprts);
  static int size(CudaMparticles* cmprts);
  static std::vector<uint> sizeByPatch(const CudaMparticles* cmprts);
  static void inject(CudaMparticles* cmprts, const std::vector<Particle>& buf,
		     const std::vector<uint>& buf_n_by_patch);

  static std::vector<uint> get_offsets(const CudaMparticles* cmprts);
  static std::vector<Particle> get_particles(const CudaMparticles* cmprts);
  static Particle get_particle(const CudaMparticles* cmprts, int p, int n);

  static bool check_after_push(CudaMparticles* cmprts);
  static void dump(const CudaMparticles* cmprts, const std::string& filename);

  static bool need_reorder(CudaMparticles* cmprts);
};

