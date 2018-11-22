
#pragma once

#include "particles.hxx"
#include "particle_cuda.hxx"
#include "particle_indexer.hxx"
#include "mparticles_patch_cuda.hxx"
#include "injector_buffered.hxx"
#include "bs.hxx"

// ======================================================================
// PatchCuda

template<typename _Mparticles>
struct ConstPatchCuda_
{
  using Mparticles = _Mparticles;
  using particle_t = typename Mparticles::particle_t;
  using real_t = typename particle_t::real_t;
  using Real3 = typename particle_t::Real3;
    
  struct const_accessor
  {
    using Double3 = Vec3<double>;
      
    const_accessor(const particle_t& prt, const ConstPatchCuda_& patch)
      : prt_{prt}, patch_{patch}
    {}

    Real3 x()   const { return prt_.x(); }
    Real3 u()   const { return prt_.u(); }
    real_t w()  const { return prt_.qni_wni() / patch_.grid().kinds[prt_.kind()].q; }
    real_t qni_wni() const { return prt_.qni_wni(); }
    int kind()  const { return prt_.kind(); }
      
    Double3 position() const
    {
      auto& patch = patch_.grid().patches[patch_.p_];
	
      return patch.xb + Double3(prt_.x());
    }
    
  private:
    particle_t prt_;
    const ConstPatchCuda_ patch_;
  };
  
  struct const_accessor_range
  {
    struct const_iterator : std::iterator<std::random_access_iterator_tag,
					  const_accessor,  // value type
					  ptrdiff_t,       // difference type
					  const_accessor*, // pointer type
					  const_accessor&> // reference type
      
    {
      const_iterator(const ConstPatchCuda_& patch, uint n)
	: patch_{patch}, n_{n}
      {}
	
      bool operator==(const_iterator other) const { return n_ == other.n_; }
      bool operator!=(const_iterator other) const { return !(*this == other); }
	
      const_iterator& operator++() { n_++; return *this; }
      const_iterator operator++(int) { auto retval = *this; ++(*this); return retval; }
      const_accessor operator*() { return {patch_.get_particle(n_), patch_}; }
	
    private:
      const ConstPatchCuda_ patch_;
      uint n_;
    };
    
    const_accessor_range(const ConstPatchCuda_& patch)
      : patch_(patch)
    {}

    const_iterator begin() const { return {patch_, 0}; }
    const_iterator end()   const { return {patch_, patch_.size()}; };
      
  private:
    const ConstPatchCuda_ patch_;
  };

  ConstPatchCuda_(Mparticles& mprts, int p)
    : mprts_(mprts), p_(p)
  {}
  
  const_accessor_range get() const { return {*this}; }

  const Grid_t& grid() const { return mprts_.grid(); }

  particle_t get_particle(int n) const
  {
    uint off = mprts_.start(p_);
    auto cprts = mprts_.get_particles(off + n, off + n + 1);
    return cprts[0];
  }

  uint size() const
  {
    uint n_prts_by_patch[grid().n_patches()];
    mprts_.get_size_all(n_prts_by_patch);
    return n_prts_by_patch[p_];
  }

protected:
  Mparticles& mprts_;
  int p_;
};

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
  
  struct Patch : ConstPatchCuda_<MparticlesCuda>
  {
    using Base = ConstPatchCuda_<MparticlesCuda>;
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
  uint start(int p) const;
  bool check_after_push();

  // ----------------------------------------------------------------------
  // facility to access particles without conversion,
  // mostly for debugging (?)

  std::vector<particle_t> get_particles(int beg, int end) const;

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

