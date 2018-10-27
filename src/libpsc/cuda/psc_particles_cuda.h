
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

// ======================================================================
// particle_cuda_t

using particle_cuda_t = psc_particle<float>;

// ======================================================================
// cuda_mparticles_prt

struct cuda_mparticles_prt
{
  using real_t = float;
  using Real3 = Vec3<real_t>;

  cuda_mparticles_prt() = default; // FIXME? needed to use std::vector

  cuda_mparticles_prt(Real3 x, Real3 p, real_t w, int kind)
    : x(x), p(p), w(w), kind(kind)
  {}

  bool operator==(const cuda_mparticles_prt& other) const
  {
    return (x == other.x && w == other.w &&
	    p == other.p && kind == other.kind);
  }

  bool operator!=(const cuda_mparticles_prt& other) const { return !(*this == other); }

  
  Real3 x;
  real_t w;
  Real3 p; 
  int kind;
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
  using particle_t = particle_cuda_t;
  using real_t = particle_t::real_t;
  using Real3 = Vec3<real_t>;
  using buf_t = std::vector<particle_cuda_t>;
  using CudaMparticles = cuda_mparticles<BS>;

  using is_cuda = std::true_type;
  
  MparticlesCuda(const Grid_t& grid);
  ~MparticlesCuda();

  int get_n_prts() const override;
  void get_size_all(uint *n_prts_by_patch) const override;
  void reserve_all(const uint *n_prts_by_patch) override;
  void resize_all(const uint *n_prts_by_patch) override;
  void reset(const Grid_t& grid) override;

  void inject_buf(const cuda_mparticles_prt *buf, const uint *buf_n_by_patch);
  void inject_buf(const particle_inject *buf, const uint *buf_n_by_patch);
  void dump(const std::string& filename);
  uint start(int p) const;
  bool check_after_push();

  // ----------------------------------------------------------------------
  // facility to access particles without conversion,
  // mostly for debugging (?)

  std::vector<cuda_mparticles_prt> get_particles(int beg, int end) const;
  std::vector<cuda_mparticles_prt> get_particles(int p) const;

  void define_species(const char *name, double q, double m,
		      double max_local_np, double max_local_nm,
		      double sort_interval, double sort_out_of_place)
  {}
  
  static const Convert convert_to_, convert_from_;
  const Convert& convert_to() override { return convert_to_; }
  const Convert& convert_from() override { return convert_from_; }

  CudaMparticles* cmprts() { return cmprts_; }

  struct patch_t
  {
    // ----------------------------------------------------------------------
    // injector
    //
    // caches injected particles for all patches before actually transferring them to
    // the GPU
    //
    // It expects that an injector is constructed and then destructed for every patch
    // in order, and the real action occurs only when the last patch instance is destructed
    
    struct injector
    {
      injector(patch_t& patch)
	: patch_{patch},
	  n_prts_{0}
      {
	auto& mprts = patch_.mp_;
	//mprintf("injector p = %d/%d\n", patch_.p_, mprts.n_patches());
	// make sure we're constructed/destructed in order, one patch at a time
	assert(patch_.p_ == mprts.injector_n_prts_by_patch.size());
      }

      ~injector()
      {
	auto& mprts = patch_.mp_;
	//mprintf("~injector p = %d/%d\n", patch_.p_, mprts.n_patches());
	mprts.injector_n_prts_by_patch.push_back(n_prts_);
	if (patch_.p_ == mprts.n_patches() - 1) {
	  mprts.inject_buf(mprts.injector_buf.data(), mprts.injector_n_prts_by_patch.data());
	  mprts.injector_n_prts_by_patch.clear();
	  mprts.injector_buf.clear();
	}
      }
      
      void operator()(const particle_inject& new_prt)
      {
	auto& mprts = patch_.mp_;
	mprts.injector_buf.push_back(new_prt);
	n_prts_++;
      }
      
    private:
      patch_t& patch_;
      uint n_prts_;
    };

    struct const_accessor
    {
    using Double3 = Vec3<double>;
      
      const_accessor(const cuda_mparticles_prt& prt, const patch_t& prts)
	: prt_{prt}, prts_{prts}
      {}
      
      Real3 u()   const { return prt_.p; }
      real_t w()  const { return prt_.w; }
      int kind()  const { return prt_.kind; }
      
      Double3 position() const
      {
	auto& patch = prts_.mp_.grid().patches[prts_.p_];
	
	return patch.xb + Double3{prt_.x};
      }
    
    private:
      cuda_mparticles_prt prt_;
      const patch_t& prts_;
    };
  
    struct const_accessor_range
    {
      struct const_iterator : std::iterator<std::random_access_iterator_tag,
					    const_accessor,  // value type
					    ptrdiff_t,       // difference type
					    const_accessor*, // pointer type
					    const_accessor&> // reference type
      
      {
	const_iterator(const patch_t& prts, uint n)
	  : prts_{prts}, n_{n}
	{}
	
	bool operator==(const_iterator other) const { return n_ == other.n_; }
	bool operator!=(const_iterator other) const { return !(*this == other); }
	
	const_iterator& operator++() { n_++; return *this; }
	const_iterator operator++(int) { auto retval = *this; ++(*this); return retval; }
	const_accessor operator*() { return {prts_.get_particle(n_), prts_}; }
	
      private:
	const patch_t& prts_;
	uint n_;
      };
    
      const_accessor_range(const patch_t& prts)
	: prts_{prts}
      {}

      const_iterator begin() const { return {prts_, 0}; }
      const_iterator end()   const { return {prts_, prts_.size()}; };
      
    private:
      const patch_t& prts_;
    };

    patch_t(MparticlesCuda& mp, int p)
      : mp_(mp), p_(p)
    {}

    const ParticleIndexer<real_t>& particleIndexer() const { return mp_.pi_; }

    cuda_mparticles_prt get_particle(int n) const
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

    injector injector() { return {*this}; }
    const_accessor_range get() const { return {*this}; }

  private:
    MparticlesCuda& mp_;
    int p_;
  };

  patch_t operator[](int p) { return patch_t{*this, p}; }

private:
  CudaMparticles* cmprts_;
  ParticleIndexer<real_t> pi_;

  std::vector<uint> injector_n_prts_by_patch;
  std::vector<particle_inject> injector_buf;
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
