
#pragma once

// FIXME?  we have two pretty similar versions of ConstPatchCuda here,
// and probably only one should survive.  one version copies from
// device "on demand", whereas the other one copies a whole patch
// worth of data just in case.  Obviously, there can be use cases for
// either, but this needs some thinking as to what's actually needed
// for this code.  The on-demand version might be useful to serve as a
// template for a modifiable accesser, if we ever want to go there.

// ======================================================================
// ConstPatchCuda

template<typename _Mparticles>
struct ConstPatchCuda
{
  using Mparticles = _Mparticles;
  using particle_t = typename Mparticles::particle_t;
  using real_t = typename particle_t::real_t;
  using Real3 = typename particle_t::Real3;

  struct const_accessor
  {
    using Double3 = Vec3<double>;
    
    const_accessor(const particle_t& prt, const Mparticles& mprts, int p)
      : prt_{prt}, mprts_{mprts}, p_{p}
    {}

    Real3 x()   const { return prt_.x(); }
    Real3 u()   const { return prt_.u(); }
    real_t w()  const { return prt_.qni_wni() / mprts_.grid().kinds[prt_.kind()].q; }
    real_t qni_wni() const { return prt_.qni_wni(); }
    int kind()  const { return prt_.kind(); }
    
    Double3 position() const
    {
      auto& patch = mprts_.grid().patches[p_];
      
      return patch.xb + Double3(prt_.x());
    }

    operator const particle_t& () const { return prt_; }
    
  private:
    const particle_t& prt_;
    const Mparticles& mprts_;
    const int p_;
  };
  
  struct const_accessor_range
  {
    struct const_iterator : std::iterator<std::random_access_iterator_tag,
      const_accessor,  // value type
      ptrdiff_t,       // difference type
      const_accessor*, // pointer type
      const_accessor&> // reference type
    
    {
      const_iterator(const const_accessor_range& range, uint n)
	: range_{range}, n_{n}
      {}
      
      bool operator==(const_iterator other) const { return n_ == other.n_; }
      bool operator!=(const_iterator other) const { return !(*this == other); }
      
      const_iterator& operator++() { n_++; return *this; }
      const_iterator operator++(int) { auto retval = *this; ++(*this); return retval; }
      const_accessor operator*() { return {range_.data_[n_], range_.mprts_, range_.p_}; }
      
    private:
      const const_accessor_range range_;
      uint n_;
    };
    
    const_accessor_range(const Mparticles& mprts, int p)
      : mprts_{mprts}, p_{p}, data_{const_cast<Mparticles&>(mprts).get_particles(p)}
    // FIXME, const hacking around reorder may change state...
    {}

    const_iterator begin() const { return {*this, 0}; }
    const_iterator end()   const { return {*this, uint(data_.size())}; }
    uint size() const { return data_.size(); }
    
  private:
    const Mparticles& mprts_;
    int p_;
    const std::vector<particle_t> data_;
  };

  ConstPatchCuda(const Mparticles& mprts, int p)
    : mprts_{mprts}, p_(p)
  {}
  
  const Grid_t& grid() const { return mprts_.grid(); }

  uint size() const
  {
    uint n_prts_by_patch[grid().n_patches()];
    mprts_.get_size_all(n_prts_by_patch);
    return n_prts_by_patch[p_];
  }

protected:
  const Mparticles& mprts_;
  int p_;
};

// ======================================================================
// ConstPatchCuda_

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
      const_accessor operator*() { return {patch_.mprts_.get_particle(patch_.p_, n_), patch_}; }
	
    private:
      const ConstPatchCuda_ patch_;
      uint n_;
    };
    
    const_accessor_range(const ConstPatchCuda_& patch)
      : patch_(patch)
    {}

    const_iterator begin() const { return {patch_, 0}; }
    const_iterator end()   const { return {patch_, patch_.size()}; };
    uint size() const { return patch_.size(); }
      
  private:
    const ConstPatchCuda_ patch_;
  };

  ConstPatchCuda_(Mparticles& mprts, int p)
    : mprts_(mprts), p_(p)
  {}
  
  const Grid_t& grid() const { return mprts_.grid(); }

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

template<typename Mparticles>
struct ConstAccessorCuda
{
  ConstAccessorCuda(Mparticles& mprts)
    : mprts_{mprts}
  {}

  typename Mparticles::Patch::const_accessor_range operator[](int p)
  {
    return {mprts_, p};
  }

private:
  Mparticles& mprts_;
};

