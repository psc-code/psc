
#ifndef PSC_PARTICLES_CUDA_H
#define PSC_PARTICLES_CUDA_H

#include "particles.hxx"
#include "particles_traits.hxx"
#include "psc_bits.h"
#include "particles_simple.hxx"

#include <vector>

struct BS144
{
  using x = std::integral_constant<unsigned int, 1>;
  using y = std::integral_constant<unsigned int, 4>;
  using z = std::integral_constant<unsigned int, 4>;
};

struct BS444
{
  using x = std::integral_constant<unsigned int, 4>;
  using y = std::integral_constant<unsigned int, 4>;
  using z = std::integral_constant<unsigned int, 4>;
};

using DParticleCuda = ParticleSimple<float>;

template<typename BS>
struct cuda_mparticles;

// ----------------------------------------------------------------------
// InjectorCuda
//
// caches injected particles for all patches before actually transferring them to
// the GPU
//
// It expects that an injector is constructed and then destructed for every patch
// in order, and the real action occurs only when the last patch instance is destructed

// ----------------------------------------------------------------------
// InjectorCuda_

template<typename Patch>
struct InjectorCuda_
{
  using particle_t = typename Patch::Mparticles::particle_t;
  using real_t = typename particle_t::real_t;
  using Real3 = typename particle_t::Real3;
  using Double3 = Vec3<double>;
  
  InjectorCuda_(const Patch& patch)
  : patch_{patch},
    n_prts_{0}
  {
    assert(patch_.p_ == patch_.cmprts_.injector_n_prts_by_patch_.size());
  }
  
  ~InjectorCuda_()
  {
    auto& cmprts = patch_.cmprts_;
    cmprts.injector_n_prts_by_patch_.push_back(n_prts_);
    if (patch_.p_ == cmprts.n_patches - 1) {
      cmprts.inject(cmprts.injector_buf_, cmprts.injector_n_prts_by_patch_);
      cmprts.injector_n_prts_by_patch_.clear();
      cmprts.injector_buf_.clear();
    }
  }
  
  void raw(const particle_t& prt)
  {
    patch_.cmprts_.injector_buf_.push_back(prt);
    n_prts_++;
  }
  
  void operator()(const particle_inject& new_prt)
  {
    auto& cmprts = patch_.cmprts_;
    auto& patch = cmprts.grid_.patches[patch_.p_];
    auto x = Double3::fromPointer(new_prt.x) - patch.xb;
    auto prt = particle_t{Real3(x), Real3(Double3::fromPointer(new_prt.u)),
			  real_t(new_prt.w), new_prt.kind};
    patch_.cmprts_.injector_buf_.push_back(prt);
    n_prts_++;
  }
  
  void operator()(const std::vector<particle_t>& buf)
  {
    auto& injector_buf = patch_.cmprts_.injector_buf_;
    injector_buf.insert(injector_buf.end(), buf.begin(), buf.end());
    n_prts_ += buf.size();
  }
  
private:
  const Patch patch_;
  uint n_prts_;
};

template<typename Patch>
struct InjectorCuda
{
  using particle_t = typename Patch::Mparticles::particle_t;
  using real_t = typename particle_t::real_t;
  using Real3 = typename particle_t::Real3;
  using Double3 = Vec3<double>;
    
  InjectorCuda(const Patch& patch)
    : patch_{patch},
      n_prts_{0}
  {
    auto& mprts = patch_.mp_;
    //mprintf("injector p = %d/%d\n", patch_.p_, mprts.n_patches());
    // make sure we're constructed/destructed in order, one patch at a time
    assert(patch_.p_ == mprts.injector_n_prts_by_patch.size());
  }
  
  ~InjectorCuda()
  {
    auto& mprts = patch_.mp_;
    //mprintf("~injector p = %d/%d\n", patch_.p_, mprts.n_patches());
    mprts.injector_n_prts_by_patch.push_back(n_prts_);
    if (patch_.p_ == mprts.n_patches() - 1) {
      mprts.inject(mprts.injector_buf, mprts.injector_n_prts_by_patch);
      mprts.injector_n_prts_by_patch.clear();
      mprts.injector_buf.clear();
    }
  }
  
  void operator()(const particle_inject& new_prt)
  {
    auto& mprts = patch_.mp_;
    auto& patch = mprts.grid().patches[patch_.p_];
    auto x = Double3::fromPointer(new_prt.x) - patch.xb;
    mprts.injector_buf.push_back({Real3(x), Real3(Double3::fromPointer(new_prt.u)),
	  real_t(new_prt.w), new_prt.kind});
    n_prts_++;
  }
  
private:
  const Patch patch_;
  uint n_prts_;
};

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
  std::vector<particle_t> get_particles(int p) const;

  void define_species(const char *name, double q, double m,
		      double max_local_np, double max_local_nm,
		      double sort_interval, double sort_out_of_place)
  {}
  
  static const Convert convert_to_, convert_from_;
  const Convert& convert_to() override { return convert_to_; }
  const Convert& convert_from() override { return convert_from_; }

  CudaMparticles* cmprts() { return cmprts_; }

  struct Patch
  {
    using Mparticles = MparticlesCuda;
    using Injector = InjectorCuda<Patch>;
    friend Injector;
    
    struct const_accessor
    {
      using Double3 = Vec3<double>;
      
      const_accessor(const particle_t& prt, const Patch& prts)
	: prt_{prt}, prts_{prts}
      {}
      
      Real3 u()   const { return prt_.u(); }
      real_t w()  const { return prt_.qni_wni() / prts_.mp_.grid().kinds[prt_.kind()].q; }
      int kind()  const { return prt_.kind(); }
      
      Double3 position() const
      {
	auto& patch = prts_.mp_.grid().patches[prts_.p_];
	
	return patch.xb + Double3(prt_.x());
      }
    
    private:
      particle_t prt_;
      const Patch prts_;
    };
  
    struct const_accessor_range
    {
      struct const_iterator : std::iterator<std::random_access_iterator_tag,
					    const_accessor,  // value type
					    ptrdiff_t,       // difference type
					    const_accessor*, // pointer type
					    const_accessor&> // reference type
      
      {
	const_iterator(const Patch& prts, uint n)
	  : prts_{prts}, n_{n}
	{}
	
	bool operator==(const_iterator other) const { return n_ == other.n_; }
	bool operator!=(const_iterator other) const { return !(*this == other); }
	
	const_iterator& operator++() { n_++; return *this; }
	const_iterator operator++(int) { auto retval = *this; ++(*this); return retval; }
	const_accessor operator*() { return {prts_.get_particle(n_), prts_}; }
	
      private:
	const Patch prts_;
	uint n_;
      };
    
      const_accessor_range(const Patch& prts)
	: prts_(prts)
      {}

      const_iterator begin() const { return {prts_, 0}; }
      const_iterator end()   const { return {prts_, prts_.size()}; };
      
    private:
      const Patch prts_;
    };

    Patch(MparticlesCuda& mp, int p)
      : mp_(mp), p_(p)
    {}

    const ParticleIndexer<real_t>& particleIndexer() const { return mp_.pi_; }

    particle_t get_particle(int n) const
    {
      uint off = mp_.start(p_);
      auto cprts = mp_.get_particles(off + n, off + n + 1);
      return cprts[0];
    }

    uint size() const
    {
      uint n_prts_by_patch[mp_.grid().n_patches()];
      mp_.get_size_all(n_prts_by_patch);
      return n_prts_by_patch[p_];
    }

    Injector injector() { return {*this}; }
    const_accessor_range get() const { return {*this}; }

  private:
    MparticlesCuda& mp_;
    int p_;
  };

  friend typename Patch::Injector;
  
  Patch operator[](int p) { return {*this, p}; }

private:
  CudaMparticles* cmprts_;
  ParticleIndexer<real_t> pi_;

  std::vector<uint> injector_n_prts_by_patch;
  std::vector<particle_t> injector_buf;
};

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


#endif
