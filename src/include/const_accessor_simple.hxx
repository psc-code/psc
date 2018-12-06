
#pragma once

// ======================================================================
// AccessorSimple

template<typename Mparticles>
struct AccessorSimple
{
  using Particle = typename Mparticles::Particle;
  using real_t = typename Mparticles::real_t;
  using Patch = typename Mparticles::Patch;
  using Real3 = Vec3<real_t>;
  
  AccessorSimple(Particle& prt, const Patch& prts)
    : prt_{prt},
      prts_{prts}
  {}
  
  real_t q() const { return prts_.prt_qni(prt_); }
  real_t m() const { return prts_.prt_mni(prt_); }
  
  real_t  u(int d) const { return prt_.u()[d]; }
  real_t& u(int d)       { return prt_.u()[d]; }
  
private:
  Particle& prt_;
  const Patch& prts_;
};

// ======================================================================
// ConstAcessorSimple

template<typename Mparticles>
struct ConstAccessorSimple
{
  using Particle = typename Mparticles::Particle;
  using MparticlesPatch = typename Mparticles::Patch;
  using Double3 = Vec3<double>;
  
  struct const_accessor
  {
    using real_t = typename Mparticles::real_t;
    using Real3 = Vec3<real_t>;
    
    const_accessor(const Particle& prt, const MparticlesPatch& prts)
      : prt_{prt}, prts_{prts}
    {}

    Real3 x()   const { return prt_.x(); }
    Real3 u()   const { return prt_.u(); }
    real_t w()  const { return prt_.qni_wni() / q(); }
    real_t qni_wni() const { return prt_.qni_wni(); }
    real_t q()  const { return prts_.grid().kinds[kind()].q; }
    real_t m()  const { return prts_.grid().kinds[kind()].m; }
    int kind()  const { return prt_.kind(); }

    Double3 position() const
    {
      auto& patch = prts_.grid().patches[prts_.p()]; // FIXME, generally, it'd be nice to have a better way to get this

      return patch.xb +	Double3(prt_.x());
    }
    
  private:
    const Particle& prt_;
    const MparticlesPatch& prts_;
  };
  
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

      const_iterator& operator++()  { n_++; return *this; }
      const_iterator operator++(int) { auto retval = *this; ++(*this); return retval; }
      const_accessor operator*() { return {prts_[n_], prts_}; }

    private:
      const MparticlesPatch& prts_;
      uint n_;
    };
    
    Patch(const MparticlesPatch& prts)
      : prts_{prts}
    {}

    const_iterator begin() const { return {prts_, 0}; }
    const_iterator end()   const { return {prts_, prts_.size()}; }
    uint size() const { return prts_.size(); }
    const_accessor operator[](int n) const { return {prts_[n], prts_}; }

  private:
    const MparticlesPatch& prts_;
  };

  ConstAccessorSimple(Mparticles& mprts)
    : mprts_{mprts}
  {}

  Patch operator[](int p) { return {mprts_[p]}; }

private:
  Mparticles& mprts_;
};

