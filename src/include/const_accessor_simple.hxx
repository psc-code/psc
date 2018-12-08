
#pragma once

// ======================================================================
// ParticleProxySimple

template<typename Mparticles>
struct ParticleProxySimple
{
  using Particle = typename Mparticles::Particle;
  using real_t = typename Mparticles::real_t;
  using Patch = typename Mparticles::Patch;
  using Real3 = Vec3<real_t>;
  
  ParticleProxySimple(Particle& prt, const Mparticles& mprts, int p)
    : prt_{prt}, mprts_{mprts}, p_{p}
  {}
  
  real_t  x(int d) const { return prt_.x()[d]; }
  real_t& x(int d)       { return prt_.x()[d]; }
  
  real_t  u(int d) const { return prt_.u()[d]; }
  real_t& u(int d)       { return prt_.u()[d]; }
  
  // FIXME, grid is always double precision, so this will switch precision
  // where not desired. should use same info stored elsewhere at right precision
  real_t w()  const { return prt_.qni_wni() / q(); }
  real_t q()  const { return mprts_.grid().kinds[kind()].q; }
  real_t m()  const { return mprts_.grid().kinds[kind()].m; }
  int kind()  const { return prt_.kind(); }

  int validCellIndex() const { return mprts_[p_].validCellIndex(prt_); }
  
  friend void swap(ParticleProxySimple<Mparticles> a, ParticleProxySimple<Mparticles> b)
  {
    using std::swap;
    swap(a.prt_, b.prt_);
  }

private:
  Particle& prt_;
  const Mparticles& mprts_;
  int p_;
};

// ======================================================================
// ConstParticleAccessorSimple

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
// AccessorPatchSimple

template<typename AccessorSimple>
struct AccessorPatchSimple
{
  using Mparticles = typename AccessorSimple::Mparticles;
  using Accessor = ParticleProxySimple<Mparticles>;
  using MparticlesPatch = typename Mparticles::Patch;
  
  AccessorPatchSimple(AccessorSimple& accessor, int p)
    : accessor_{accessor}, p_{p}
  {}

  Accessor operator[](int n) { return {accessor_.data(p_)[n], accessor_.mprts(), p_}; }
  uint size() const { return accessor_.size(p_); }
  const Grid_t& grid() const { return accessor_.grid(); }

private:
  AccessorSimple& accessor_;
  int p_;
};

// ======================================================================
// ConstAcessorPatchSimple

template<typename ConstAccessorSimple>
struct ConstAccessorPatchSimple
{
  using Mparticles = typename ConstAccessorSimple::Mparticles;
  using accessor = ConstParticleAccessorSimple<Mparticles>;

  struct const_iterator : std::iterator<std::forward_iterator_tag,
					accessor,        // value type
					ptrdiff_t,       // difference type
					const accessor*, // pointer type
					const accessor&> // reference type
  
  {
    const_iterator(const ConstAccessorPatchSimple& patch, uint n)
      : patch_{patch}, n_{n}
    {}
    
    bool operator==(const_iterator other) const { return n_ == other.n_; }
    bool operator!=(const_iterator other) const { return !(*this == other); }
    
    const_iterator& operator++() { n_++; return *this; }
    const_iterator operator++(int) { auto retval = *this; ++(*this); return retval; }
    const accessor operator*() { return patch_[n_]; }
    
  private:
    const ConstAccessorPatchSimple patch_;
    uint n_;
  };
  
  ConstAccessorPatchSimple(const ConstAccessorSimple& accessor, int p)
    : accessor_{accessor}, p_{p}
  {}
  
  const_iterator begin() const { return {*this, 0}; }
  const_iterator end()   const { return {*this, size()}; }
  const accessor operator[](int n) const { return {accessor_.data(p_)[n], accessor_.mprts(), p_}; }
  uint size() const { return accessor_.size(p_); }
  
private:
  const ConstAccessorSimple& accessor_;
  const int p_;
};

// ======================================================================
// ConstAcessorSimple

template<typename _Mparticles>
struct ConstAccessorSimple
{
  using Mparticles = _Mparticles;
  using Patch = ConstAccessorPatchSimple<ConstAccessorSimple>;
  
  ConstAccessorSimple(Mparticles& mprts)
    : mprts_{mprts}
  {}

  Patch operator[](int p) const { return {*this, p}; }
  const Mparticles& mprts() const { return mprts_; }
  uint size(int p) const { return mprts_[p].size(); }
  typename Mparticles::Patch::iterator data(int p) const { return mprts_[p].begin(); }
  const Grid_t& grid() const { return mprts_.grid(); }

private:
  Mparticles& mprts_;
};

// ======================================================================
// AcessorSimple

template<typename _Mparticles>
struct AccessorSimple
{
  using Mparticles = _Mparticles;
  using Patch = AccessorPatchSimple<AccessorSimple>;
  
  AccessorSimple(Mparticles& mprts)
    : mprts_{mprts}
  {}

  Patch operator[](int p) { return {*this, p}; }
  const Mparticles& mprts() const { return mprts_; }
  Mparticles& mprts() { return mprts_; }
  uint size(int p) const { return mprts_[p].size(); }
  typename Mparticles::Patch::iterator data(int p) { return mprts_[p].begin(); }
  const Grid_t& grid() const { return mprts_.grid(); }

private:
  Mparticles& mprts_;
};

