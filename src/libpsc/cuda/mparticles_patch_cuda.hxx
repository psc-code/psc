
#pragma once

#include "const_accessor_simple.hxx"

// FIXME?  we have two pretty similar versions of ConstAccessorCuda here,
// and probably only one should survive.  one version copies from
// device "on demand", whereas the other one copies a whole patch
// worth of data just in case.  Obviously, there can be use cases for
// either, but this needs some thinking as to what's actually needed
// for this code.  The on-demand version might be useful to serve as a
// template for a modifiable accesser, if we ever want to go there.

// ======================================================================
// ConstAccessorCuda

template<typename _Mparticles>
struct ConstAccessorCuda
{
  using Mparticles = _Mparticles;
  using _Particle = typename Mparticles::Particle;
  using Patch = ConstAccessorPatchSimple<ConstAccessorCuda>;
  using Particle = typename Patch::ConstParticleProxy;
  
  ConstAccessorCuda(Mparticles& mprts)
    : mprts_{mprts}, data_{const_cast<Mparticles&>(mprts).get_particles()}, off_{mprts.get_offsets()}
  {}

  Patch operator[](int p) const 
    {
        printf("REORDER: %d\n", mprts_.need_reorder());
        if(mprts_.need_reorder())
        {
            printf("Warning: Calling accessor with unordered particles! Expect invalid results\n");
            abort();
        } 

         return {*this, p}; 
    }
  Mparticles& mprts() const { return mprts_; }
  const _Particle* data(int p) const { return &data_[off_[p]]; }
  uint size(int p) const { return off_[p+1] - off_[p]; }

private:
  Mparticles& mprts_;
  const std::vector<_Particle> data_;
  const std::vector<uint> off_;
};

// ======================================================================
// ConstParticleAccessorCuda_

template<typename Mparticles>
struct ConstParticleAccessorPatchCuda
{
  ConstParticleAccessorPatchCuda(const Mparticles& mprts, int p)
    : mprts{mprts}, p{p}
  {}

  const Grid_t& grid() const { return mprts.grid(); }

  const Mparticles& mprts;
  int p;
};

template<typename Mparticles>
struct ConstParticleAccessorCuda_
{
  using Particle = typename Mparticles::Particle;
  using real_t = typename Particle::real_t;
  using Real3 = Vec3<real_t>;
  using Double3 = Vec3<double>;
  using MparticlesPatch = ConstParticleAccessorPatchCuda<Mparticles>;

  ConstParticleAccessorCuda_(const Particle& prt, const MparticlesPatch& prts)
    : prt_{prt}, prts_{prts}
  {}
  
  Real3 x()   const { return prt_.x(); }
  Real3 u()   const { return prt_.u(); }
  real_t w()  const { return prt_.qni_wni() / prts_.grid().kinds[prt_.kind()].q; }
  real_t qni_wni() const { return prt_.qni_wni(); }
  int kind()  const { return prt_.kind(); }
  
  Double3 position() const
  {
    auto& patch = prts_.grid().patches[prts_.p];
    
    return patch.xb + Double3(prt_.x());
  }
  
  operator const Particle& () const { return prt_; }
  
private:
  Particle prt_;
  MparticlesPatch prts_;
};
  

// ======================================================================
// ConstAccessorCuda_

template<typename Mparticles>
struct ConstAccessorCuda_
{
  using Particle = typename Mparticles::Particle;
  using real_t = typename Particle::real_t;
  using Real3 = typename Particle::Real3;
  using MparticlesPatch = ConstParticleAccessorPatchCuda<Mparticles>;

  using const_accessor = ConstParticleAccessorCuda_<Mparticles>;
  
  struct Patch
  {
    struct const_iterator : std::iterator<std::random_access_iterator_tag,
					  const_accessor,  // value type
					  ptrdiff_t,       // difference type
					  const_accessor*, // pointer type
					  const_accessor&> // reference type
      
    {
      const_iterator(const MparticlesPatch& prts, uint n)
	: prts_{prts}, n_{n}
      {}
	
      bool operator==(const_iterator other) const { return n_ == other.n_; }
      bool operator!=(const_iterator other) const { return !(*this == other); }
	
      const_iterator& operator++() { n_++; return *this; }
      const_iterator operator++(int) { auto retval = *this; ++(*this); return retval; }
      const_accessor operator*() { return {const_cast<Mparticles&>(prts_.mprts).get_particle(prts_.p, n_), prts_}; } // FIXME constness
	
    private:
      const MparticlesPatch prts_;
      uint n_;
    };
    
    Patch(const Mparticles& mprts, int p)
      : prts_{mprts, p}
    {}

    const_iterator begin() const { return {prts_, 0}; }
    const_iterator end()   const { return {prts_, size()}; };
    const_accessor operator[](int n) const {return {const_cast<Mparticles&>(prts_.mprts).get_particle(prts_.p, n), prts_}; }

    uint size() const
    {
      auto n_prts_by_patch = const_cast<Mparticles&>(prts_.mprts).sizeByPatch();
      return n_prts_by_patch[prts_.p];
    }

      
  private:
    const MparticlesPatch prts_;
  };

  ConstAccessorCuda_(Mparticles& mprts)
    : mprts_{mprts}
  {}

  Patch operator[](int p) { return {mprts_, p}; }

private:
  Mparticles& mprts_;
};

// ======================================================================
// ConstPatchCuda

template<typename _Mparticles>
struct ConstPatchCuda
{
  using Mparticles = _Mparticles;

  ConstPatchCuda(const Mparticles& mprts, int p)
    : mprts_{mprts}, p_(p)
  {}
  
private:
  const Mparticles& mprts_;
  int p_;
};

