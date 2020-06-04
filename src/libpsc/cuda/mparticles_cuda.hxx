
#pragma once

#include "particles.hxx"
#include "particle_cuda.hxx"
#include "particle_indexer.hxx"
#include "mparticles_patch_cuda.hxx"
#include "injector_buffered.hxx"
#include "cuda_mparticles_iface.hxx"
#include "bs.hxx"
#include "UniqueIdGenerator.h"

// ======================================================================
// MparticlesCuda

template<typename _BS>
struct MparticlesCuda : MparticlesBase
{
  using Self = MparticlesCuda;
  using BS = _BS;
  using real_t = float;
  using Particle = ParticleSimple<real_t>;
  using Real3 = Vec3<real_t>;
  using CudaMparticles = cuda_mparticles<BS>;
  using BndpParticle = DParticleCuda;
  using BndBuffer = std::vector<BndpParticle>;
  using BndBuffers = std::vector<BndBuffer>;
  using ConstAccessor = ConstAccessorCuda<MparticlesCuda>;
  using is_cuda = std::true_type;
  
  using Patch = ConstPatchCuda<MparticlesCuda>;
  using Iface = cuda_mparticles_iface<BS>;

  MparticlesCuda(const Grid_t& grid)
    : MparticlesBase(grid),
      pi_(grid),
      uid_gen(grid.comm())
  {
    cmprts_ = Iface::new_(grid);
  }
  
  ~MparticlesCuda() { Iface::delete_(cmprts_); }

  void reset(const Grid_t& grid) override
  {
    this->~MparticlesCuda();
    new(this) MparticlesCuda(grid);
  }

  int size() const override { return Iface::size(cmprts_); }
  std::vector<uint> sizeByPatch() const override { return Iface::sizeByPatch(cmprts_); }

  void inject(const std::vector<Particle>& buf, const std::vector<uint>& buf_n_by_patch) { Iface::inject(cmprts_, buf, buf_n_by_patch); }
  void dump(const std::string& filename) { Iface::dump(cmprts_, filename); }
  bool check_after_push() { return Iface::check_after_push(cmprts_); }
  void clear() { Iface::clear(cmprts_); }

  std::vector<uint> get_offsets() const { return Iface::get_offsets(cmprts_); }
  std::vector<Particle> get_particles() const { return Iface::get_particles(cmprts_); }
  Particle get_particle(int p, int n) const { return Iface::get_particle(cmprts_, p, n); }
  const ParticleIndexer<real_t>& particleIndexer() const { return pi_; }
  void define_species(const char *name, double q, double m,
		      double max_local_np, double max_local_nm,
		      double sort_interval, double sort_out_of_place)
  {}
  
  static const Convert convert_to_, convert_from_;
  const Convert& convert_to() override { return convert_to_; }
  const Convert& convert_from() override { return convert_from_; }

  CudaMparticles* cmprts() { return cmprts_; }

  InjectorBuffered<MparticlesCuda> injector() { return {*this}; }
  ConstAccessorCuda<MparticlesCuda> accessor() const { return {const_cast<MparticlesCuda&>(*this)}; } // FIXME cast

private:
  CudaMparticles* cmprts_;
  ParticleIndexer<real_t> pi_;
public:
  psc::particle::UniqueIdGenerator uid_gen;
};

