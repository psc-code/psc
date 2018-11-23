
#pragma once

#include "particles.hxx"
#include "particle_cuda.hxx"
#include "particle_indexer.hxx"
#include "mparticles_patch_cuda.hxx"
#include "injector_buffered.hxx"
#include "bs.hxx"

template<typename BS>
struct cuda_mparticles;

// ======================================================================
// MparticlesCuda

template<typename _BS>
struct MparticlesCuda : MparticlesBase
{
  using Self = MparticlesCuda;
  using BS = _BS;
  using real_t = float;
  using particle_t = ParticleSimple<real_t>;
  using Real3 = Vec3<real_t>;
  using BndpParticle = DParticleCuda;
  using buf_t = std::vector<BndpParticle>;
  using CudaMparticles = cuda_mparticles<BS>;

  using is_cuda = std::true_type;
  
  struct Patch : ConstPatchCuda<MparticlesCuda>
  {
    using Base = ConstPatchCuda<MparticlesCuda>;
    using Base::Base;
    using Base::mprts_;
    using Base::p_;
    using Base::grid;

    const ParticleIndexer<real_t>& particleIndexer() const { return mprts_.pi_; }
  };

  MparticlesCuda(const Grid_t& grid);
  ~MparticlesCuda();

  int get_n_prts() const override;
  void get_size_all(uint *n_prts_by_patch) const override;
  void reset(const Grid_t& grid) override;

  void inject(const std::vector<particle_t>& buf, const std::vector<uint>& buf_n_by_patch);
  void dump(const std::string& filename);
  bool check_after_push();

  std::vector<particle_t> get_particles(int p) const;
  particle_t get_particle(int p, int n) const;
  
  void define_species(const char *name, double q, double m,
		      double max_local_np, double max_local_nm,
		      double sort_interval, double sort_out_of_place)
  {}
  
  static const Convert convert_to_, convert_from_;
  const Convert& convert_to() override { return convert_to_; }
  const Convert& convert_from() override { return convert_from_; }

  CudaMparticles* cmprts() { return cmprts_; }

  Patch operator[](int p) { return {*this, p}; }
  InjectorBuffered<MparticlesCuda> injector() { return {*this}; }

private:
  CudaMparticles* cmprts_;
  ParticleIndexer<real_t> pi_;
};

// FIXME
template<>
struct Mparticles_traits<MparticlesCuda<BS144>>
{
  static constexpr const char* name = "cuda";
  static MPI_Datatype mpi_dtype() { return MPI_FLOAT; }
};


template<>
struct Mparticles_traits<MparticlesCuda<BS444>>
{
  static constexpr const char* name = "cuda444";
  static MPI_Datatype mpi_dtype() { return MPI_FLOAT; }
};

