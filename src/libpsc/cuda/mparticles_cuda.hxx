
#pragma once

#include "particles.hxx"
#include "particle_cuda.hxx"
#include "particle_indexer.hxx"
#include "mparticles_patch_cuda.hxx"
#include "injector_buffered.hxx"
#include "bs.hxx"
#include "UniqueIdGenerator.h"

#include "cuda_mparticles.hxx"

// ======================================================================
// MparticlesCuda

template <typename _BS>
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

  MparticlesCuda(const Grid_t& grid)
    : MparticlesBase(grid), pi_(grid), uid_gen(grid.comm())
  {
    cmprts_ = new CudaMparticles(grid);
  }

  ~MparticlesCuda() { delete cmprts_; }

  void reset(const Grid_t& grid) override
  {
    this->~MparticlesCuda();
    new (this) MparticlesCuda(grid);
  }

  int size() const override { return cmprts_->size(); }

  std::vector<uint> sizeByPatch() const override
  {
    return cmprts_->sizeByPatch();
  }

  void inject(const std::vector<Particle>& buf,
              const std::vector<uint>& buf_n_by_patch)
  {
    cmprts_->inject(buf, buf_n_by_patch);
  }

  double mem_fraction() const { return cmprts_->mem_fraction(); }

  void dump(const std::string& filename) { cmprts_->dump(filename); }

  bool check_after_push() { return cmprts_->check_bidx_after_push(); }

  void clear() { cmprts_->clear(); }

  std::vector<uint> get_offsets() const { return cmprts_->get_offsets(); }

  std::vector<Particle> get_particles() const
  {
    return const_cast<CudaMparticles*>(cmprts_)->get_particles(); // FIXME cast
  }

  const ParticleIndexer<real_t>& particleIndexer() const { return pi_; }
  void define_species(const char* name, double q, double m, double max_local_np,
                      double max_local_nm, double sort_interval,
                      double sort_out_of_place)
  {}

  static const Convert convert_to_, convert_from_;
  const Convert& convert_to() override { return convert_to_; }
  const Convert& convert_from() override { return convert_from_; }

  CudaMparticles* cmprts() { return cmprts_; }

  InjectorBuffered<MparticlesCuda> injector() { return {*this}; }
  ConstAccessorCuda<MparticlesCuda> accessor() const
  {
    return {const_cast<MparticlesCuda&>(*this)};
  } // FIXME cast

private:
  CudaMparticles* cmprts_;
  ParticleIndexer<real_t> pi_;

public:
  psc::particle::UniqueIdGenerator uid_gen;
};
