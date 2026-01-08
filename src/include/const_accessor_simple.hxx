
#pragma once

// ======================================================================
// ParticleProxySimple

template <typename Mparticles>
struct ParticleProxySimple
{
  using Particle = typename Mparticles::Particle;
  using real_t = typename Mparticles::real_t;
  using Patch = typename Mparticles::Patch;
  using Real3 = Vec3<real_t>;

  ParticleProxySimple(Particle& prt, const Mparticles& mprts)
    : prt_{prt}, mprts_{mprts}
  {}

  Real3 x() const { return prt_.x; }
  Real3& x() { return prt_.x; }

  Real3 u() const { return prt_.u; }
  Real3& u() { return prt_.u; }

  // FIXME, grid is always double precision, so this will switch precision
  // where not desired. should use same info stored elsewhere at right precision
  real_t qni_wni() const { return prt_.qni_wni; }
  real_t w() const { return prt_.qni_wni / q(); }
  real_t q() const { return mprts_.prt_q(prt_); }
  real_t m() const { return mprts_.prt_m(prt_); }
  int kind() const { return prt_.kind; }

  friend void swap(ParticleProxySimple<Mparticles> a,
                   ParticleProxySimple<Mparticles> b)
  {
    using std::swap;
    swap(a.prt_, b.prt_);
  }

private:
  Particle& prt_;
  const Mparticles& mprts_;
};

// ======================================================================
// ConstParticleProxySimple

template <typename Mparticles>
struct ConstParticleProxySimple
{
  using Particle = typename Mparticles::Particle;
  using real_t = typename Mparticles::real_t;
  using Real3 = Vec3<real_t>;
  using Double3 = Vec3<double>;

  ConstParticleProxySimple(const Particle& prt, const Mparticles& mprts, int p)
    : prt_{prt}, mprts_{mprts}, p_{p}
  {}

  Real3 x() const { return prt_.x; }
  Real3 u() const { return prt_.u; }
  real_t w() const { return prt_.qni_wni / q(); }
  real_t qni_wni() const { return prt_.qni_wni; }
  real_t q() const { return mprts_.grid().kinds[kind()].q; }
  real_t m() const { return mprts_.grid().kinds[kind()].m; }
  int kind() const { return prt_.kind; }
  psc::particle::Id id() const { return prt_.id(); }
  psc::particle::Tag tag() const { return prt_.tag(); }

  Double3 position() const
  {
    auto& patch = mprts_.grid().patches[p_]; // FIXME, generally, it'd be nice
                                             // to have a better way to get this

    return patch.xb + Double3(prt_.x);
  }

private:
  const Particle& prt_;
  const Mparticles& mprts_;
  const int p_;
};

// ======================================================================
// AccessorPatchSimple

template <typename AccessorSimple>
struct AccessorPatchSimple
{
  using Mparticles = typename AccessorSimple::Mparticles;
  using ParticleProxy = ParticleProxySimple<Mparticles>;

  struct iterator
  {
    using iterator_category = std::forward_iterator_tag;
    using value_type = ParticleProxy;
    using difference_type = std::ptrdiff_t;
    using pointer = ParticleProxy*;
    using reference = ParticleProxy&;

    iterator(AccessorPatchSimple& patch, uint n) : patch_{patch}, n_{n} {}

    bool operator==(iterator other) const { return n_ == other.n_; }
    bool operator!=(iterator other) const { return !(*this == other); }

    iterator& operator++()
    {
      n_++;
      return *this;
    }
    iterator operator++(int)
    {
      auto retval = *this;
      ++(*this);
      return retval;
    }
    ParticleProxy operator*() { return patch_[n_]; }

  private:
    AccessorPatchSimple patch_;
    uint n_;
  };

  AccessorPatchSimple(AccessorSimple& accessor, int p)
    : accessor_{accessor}, p_{p}
  {}

  iterator begin() { return {*this, 0}; }
  iterator end() { return {*this, size()}; }
  ParticleProxy operator[](int n)
  {
    return {accessor_.data(p_)[n], accessor_.mprts()};
  }
  uint size() const { return accessor_.size(p_); }
  const Grid_t& grid() const { return accessor_.grid(); }

private:
  AccessorSimple& accessor_;
  int p_;
};

// ======================================================================
// ConstAccessorPatchSimple

template <typename ConstAccessorSimple>
struct ConstAccessorPatchSimple
{
  using Mparticles = typename ConstAccessorSimple::Mparticles;
  using ConstParticleProxy = ConstParticleProxySimple<Mparticles>;

  struct const_iterator
  {
    using iterator_category = std::forward_iterator_tag;
    using value_type = ConstParticleProxy;
    using difference_type = std::ptrdiff_t;
    using pointer = ConstParticleProxy*;
    using reference = ConstParticleProxy&;

    const_iterator(const ConstAccessorPatchSimple& patch, uint n)
      : patch_{patch}, n_{n}
    {}

    bool operator==(const_iterator other) const { return n_ == other.n_; }
    bool operator!=(const_iterator other) const { return !(*this == other); }

    const_iterator& operator++()
    {
      n_++;
      return *this;
    }
    const_iterator operator++(int)
    {
      auto retval = *this;
      ++(*this);
      return retval;
    }
    ConstParticleProxy operator*() { return patch_[n_]; }

  private:
    const ConstAccessorPatchSimple patch_;
    uint n_;
  };

  ConstAccessorPatchSimple(const ConstAccessorSimple& accessor, int p)
    : accessor_{accessor}, p_{p}
  {}

  const_iterator begin() const { return {*this, 0}; }
  const_iterator end() const { return {*this, size()}; }
  ConstParticleProxy operator[](int n) const
  {
    return {accessor_.data(p_)[n], accessor_.mprts(), p_};
  }
  uint size() const { return accessor_.size(p_); }

private:
  const ConstAccessorSimple& accessor_;
  const int p_;
};

// ======================================================================
// ConstAccessorSimple

template <typename _Mparticles>
struct ConstAccessorSimple
{
  using Mparticles = _Mparticles;
  using Patch = ConstAccessorPatchSimple<ConstAccessorSimple>;
  using Particle = typename Patch::ConstParticleProxy;

  ConstAccessorSimple(Mparticles& mprts) : mprts_{mprts} {}

  Patch operator[](int p) const { return {*this, p}; }
  const Mparticles& mprts() const { return mprts_; }
  uint size(int p) const { return mprts_.size(p); }
  typename Mparticles::Patch::iterator data(int p) const
  {
    return mprts_[p].begin();
  }
  const Grid_t& grid() const { return mprts_.grid(); }

private:
  Mparticles& mprts_;
};

// ======================================================================
// AccessorSimple

template <typename _Mparticles>
struct AccessorSimple
{
  using Mparticles = _Mparticles;
  using Patch = AccessorPatchSimple<AccessorSimple>;

  AccessorSimple(Mparticles& mprts) : mprts_{mprts} {}

  Patch operator[](int p) { return {*this, p}; }
  const Mparticles& mprts() const { return mprts_; }
  Mparticles& mprts() { return mprts_; }
  uint size(int p) const { return mprts_.size(p); }
  typename Mparticles::Patch::iterator data(int p) { return mprts_[p].begin(); }
  const Grid_t& grid() const { return mprts_.grid(); }

private:
  Mparticles& mprts_;
};
