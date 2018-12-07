
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
// ConstParticleAcessorSimple

template<typename Mparticles>
struct ConstParticleAccessorSimple
{
  using Particle = typename Mparticles::Particle;
  using real_t = typename Mparticles::real_t;
  using Real3 = Vec3<real_t>;
  using Double3 = Vec3<double>;
  
  ConstParticleAccessorSimple(const Particle& prt, const Mparticles& mprts, int p)
    : prt_{prt}, mprts_{mprts}, p_{p}
  {}
  
  Real3 x()   const { return prt_.x(); }
  Real3 u()   const { return prt_.u(); }
  real_t w()  const { return prt_.qni_wni() / q(); }
  real_t qni_wni() const { return prt_.qni_wni(); }
  real_t q()  const { return mprts_.grid().kinds[kind()].q; }
  real_t m()  const { return mprts_.grid().kinds[kind()].m; }
  int kind()  const { return prt_.kind(); }
  
  Double3 position() const
  {
    auto& patch = mprts_.grid().patches[p_]; // FIXME, generally, it'd be nice to have a better way to get this
    
    return patch.xb + Double3(prt_.x());
  }
  
private:
  const Particle& prt_;
  const Mparticles& mprts_;
  const int p_;
};
  
// ======================================================================
// ConstAcessorPatchSimple

template<typename Mparticles>
struct ConstAccessorSimple;

template<typename Mparticles>
struct ConstAccessorPatchSimple
{
  using const_accessor = ConstParticleAccessorSimple<Mparticles>;
  using ConstAccessorSimple = ConstAccessorSimple<Mparticles>;
  using MparticlesPatch = typename Mparticles::Patch;

  struct const_iterator : std::iterator<std::random_access_iterator_tag,
					const_accessor,  // value type
					ptrdiff_t,       // difference type
					const_accessor*, // pointer type
					const_accessor&> // reference type
  
  {
    const_iterator(const ConstAccessorPatchSimple& patch, uint n)
      : patch_{patch}, n_{n}
    {}
    
    bool operator==(const_iterator other) const { return n_ == other.n_; }
    bool operator!=(const_iterator other) const { return !(*this == other); }
    
    const_iterator& operator++() { n_++; return *this; }
    const_iterator operator++(int) { auto retval = *this; ++(*this); return retval; }
    const_accessor operator*() { return patch_[n_]; }
    
  private:
    const ConstAccessorPatchSimple patch_;
    uint n_;
  };
  
  ConstAccessorPatchSimple(const ConstAccessorSimple& accessor, int p)
    : prts_{accessor.mprts()[p]}
  {}
  
  const_iterator begin() const { return {*this, 0}; }
  const_iterator end()   const { return {*this, size()}; }
  const_accessor operator[](int n) const { return {prts_[n], prts_.mprts(), prts_.p()}; }
  uint size() const { return prts_.size(); }
  
private:
  const MparticlesPatch& prts_;
};

// ======================================================================
// ConstAcessorSimple

template<typename Mparticles>
struct ConstAccessorSimple
{
  using Patch = ConstAccessorPatchSimple<Mparticles>;
  
  ConstAccessorSimple(Mparticles& mprts)
    : mprts_{mprts}
  {}

  Patch operator[](int p) const { return {*this, p}; }

  const Mparticles& mprts() const { return mprts_; }

private:
  Mparticles& mprts_;
};

