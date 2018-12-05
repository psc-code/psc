
#pragma once

#include "psc_vpic_bits.h"

struct vpic_mparticles_prt
{
  float dx[3];
  int i;
  float ux[3];
  float w;
  int kind;
};

// ======================================================================
// InjectorVpic

template<typename Mparticles>
struct InjectorVpic
{
  struct Patch
  {
    Patch(Mparticles& mprts)
      : mprts_{mprts}
    {}

    void operator()(const particle_inject& prt)
    {
      const auto& vgrid = mprts_.vgrid();
      float dVi = 1.f / (vgrid.dx * vgrid.dy * vgrid.dz);
      particle_inject prt_reweighted = prt;
      prt_reweighted.w *= dVi;
      reweight(prt_reweighted);
    }
    
    void reweight(const particle_inject& prt)
    {
      mprts_.inject_particle_reweight(prt);
    }

  private:
    Mparticles& mprts_;
  };
  
  InjectorVpic(Mparticles& mprts)
    : mprts_{mprts}
  {}
  
  Patch operator[](int p) const { return {mprts_}; }

private:
  Mparticles& mprts_;
};

// ======================================================================
// ConstAccessorVpic

template<typename Mparticles>
struct ConstAccessorVpic
{
  using Particle = typename Mparticles::Particle;
  using const_iterator = typename Mparticles::const_iterator;
  using real_t = float;
  using Real3 = Vec3<real_t>;
  using Double3 = Vec3<double>;

  struct const_accessor
  {
    const_accessor(const_iterator sp, uint n)
      : sp_{sp}, n_{n}
    {}
      
    Real3 u()  const { return {prt().ux, prt().uy, prt().uz}; }
    real_t w() const { return prt().w * sp_->vgrid().dV; }
    real_t qni_wni() const { return w() * sp_->q; }
    int kind() const { return sp_->id; }
      
    Real3 x() const { return Real3(x_double()); }
      
    Double3 x_double()  const
    {
      const auto& vgrid = sp_->vgrid();
      double x0 = vgrid.x0, y0 = vgrid.y0, z0 = vgrid.z0;
      double x1 = vgrid.x1, y1 = vgrid.y1, z1 = vgrid.z1;
      double nx = vgrid.nx, ny = vgrid.ny, nz = vgrid.nz;
	
      int i = prt().i;
      int iz = i / ((nx+2) * (ny+2));
      i -= iz * ((nx+2) * (ny+2));
      int iy = i / (nx+2);
      i -= iy * (nx + 2);
      int ix = i;
	
      // adjust to 0-based (no ghost)
      ix--; iy--; iz--;
	
      // back to physical coords
      Double3 x = { ix + .5*(prt().dx+1.),
		    iy + .5*(prt().dy+1.),
		    iz + .5*(prt().dz+1.) };
      x *= (Double3{x1, y1, z1} - Double3{x0, y0, z0}) / Double3{nx, ny, nz};
	
      return x;
    }
      
    Double3 position() const
    {
      const auto& vgrid = sp_->vgrid();
      double x0 = vgrid.x0, y0 = vgrid.y0, z0 = vgrid.z0;
	
      return Double3(x_double()) + Double3{x0, y0, z0};
    }
      
  private:
    const Particle& prt() const { return sp_->p[n_]; }
      
    const_iterator sp_;
    uint n_;
  };
    
  struct Patch
  {
    struct const_iterator : std::iterator<std::random_access_iterator_tag,
					  const_accessor,  // value type
					  ptrdiff_t,       // difference type
					  const_accessor*, // pointer type
					  const_accessor&> // reference type
      
    {
      const_iterator(const Mparticles& mprts, typename Mparticles::const_iterator sp, uint n)
	: mprts_{mprts}, sp_{sp}, n_{n}
      {}
	
      bool operator==(const_iterator other) const { return sp_ == other.sp_ && n_ == other.n_; }
      bool operator!=(const_iterator other) const { return !(*this == other); }
	
      const_iterator& operator++()
      {
	n_++;
	if (n_ == sp_->np) {
	  n_ = 0;
	  ++sp_;
	}
	return *this;
      }
      
      const_iterator operator++(int) { auto retval = *this; ++(*this); return retval; }
      const_accessor operator*() { return {sp_, n_}; }
	
    private:
      const Mparticles& mprts_;
      typename Mparticles::const_iterator sp_;
      uint n_;
    };
      
    Patch(const Mparticles& mprts)
      : mprts_{mprts}
    {}
      
    const_iterator begin() const { return {mprts_, mprts_.begin(), 0}; }
    const_iterator end()   const { return {mprts_, mprts_.end(), 0}; }
    uint size() const { return mprts_.get_n_prts(); }
      
  private:
    const Mparticles& mprts_;
  };
  
  ConstAccessorVpic(Mparticles& mprts)
    : mprts_{mprts}
  {}

  Patch operator[](int p) { return {mprts_}; }

private:
  Mparticles& mprts_;
};

// ======================================================================
// MparticlesVpic_

template<typename _Particles>
struct MparticlesVpic_ : MparticlesBase, _Particles
{
  using Particles = _Particles;
  using Species = typename Particles::Species;
  using Grid = typename Particles::Grid;
  using Particle = typename Particles::Particle;
  using const_iterator = typename Particles::const_iterator;
  using real_t = float;

  // FIXME, should go away or actually made uesful
  struct Patch {};

  // ----------------------------------------------------------------------
  // ctor

  MparticlesVpic_(const Grid_t& grid, Grid* vgrid)
    : MparticlesBase(grid),
      vgrid_(vgrid)
  {
    assert(grid.n_patches() == 1);
  }

  // ----------------------------------------------------------------------
  // get_n_prts

  int get_n_prts() const override
  {
    int n_prts = 0;
    for (auto& sp : *this) {
      n_prts += sp.np;
    }
    
    return n_prts;
  }

  // ----------------------------------------------------------------------
  // get_size_all
  
  std::vector<uint> get_size_all() const override
  {
    return {uint(get_n_prts())};
  }

  // ----------------------------------------------------------------------
  // reserve_all
  //
  // This is a bit iffy, since we don't really want to reallocate stuff here,
  // at least for now, and we wouldn't be able to know how to split this into
  // the different species, anyway.

  void reserve_all(const std::vector<uint>& n_prts_by_patch)
  {
    for (int p = 0; p < n_patches(); p++) {
      int n_prts = 0, n_prts_alloced = 0;
      for (auto& sp : *this) {
	n_prts += sp.np;
	n_prts_alloced += sp.max_np;
      }
#if 0
      if (n_prts_by_patch[p] != n_prts) {
	mprintf("vpic_mparticles_reserve_all: %d (currently %d max %d)\n",
		n_prts_by_patch[p], n_prts, n_prts_alloced);
      }
#endif
      assert(n_prts_by_patch[p] <= n_prts_alloced);
    }
  }

  // ----------------------------------------------------------------------
  // resize_all
  //
  // Even more iffy, since can't really resize the per-species arrays, since we don't
  // know how the total # of particles we're given should be divided up
  
  void resize_all(const std::vector<uint>& n_prts_by_patch)
  {
    // we can't resize to the numbers given, unless it's "resize to 0", we'll just do nothing
    // The mparticles conversion function should call resize_all() itself first, resizing to
    // 0, and then using push_back, which will increase the count back to the right value

    if (n_prts_by_patch[0] == 0) {
      for (auto& sp : *this) {
	sp.np = 0;
      }
    } else {
#if 0
      int cur_n_prts_by_patch[n_patches];
      get_size_all(n_patches, cur_n_prts_by_patch);
      
      mprintf("vpic_mparticles_resize_all: ignoring %d -> %d\n",
	      cur_n_prts_by_patch[0], n_prts_by_patch[0]);
#endif
    }
  }

  void reset()
  {
    for (auto& sp : *this) {
      sp.np = 0;
    }
  }

  void push_back(const vpic_mparticles_prt *prt)
  {
    for (auto& sp : *this) {
      if (sp.id == prt->kind) {
	assert(sp.np < sp.max_np);
	// the below is inject_particle_raw()
	auto * RESTRICT p = sp.p + (sp.np++);
	p->dx = prt->dx[0]; p->dy = prt->dx[1]; p->dz = prt->dx[2]; p->i = prt->i;
	p->ux = prt->ux[0]; p->uy = prt->ux[1]; p->uz = prt->ux[2]; p->w = prt->w;
	return;
      }
    }
    mprintf("prt->kind %d not found in species list!\n", prt->kind);
    assert(0);
  }
  
  InjectorVpic<MparticlesVpic_> injector() { return {*this}; }

  ConstAccessorVpic<MparticlesVpic_> accessor() { return {*this}; }

  Species* define_species(const char *name, double q, double m,
			  double max_local_np, double max_local_nm,
			  double sort_interval, double sort_out_of_place)
  {
    // Compute a reasonble number of movers if user did not specify
    // Based on the twice the number of particles expected to hit the boundary
    // of a wpdt=0.2 / dx=lambda species in a 3x3x3 domain
    if (max_local_nm < 0) {
      max_local_nm = 2 * max_local_np / 25;
#if 0
      // FIXME, don't know MAX_PIPELINE, and that's mostly gone
      // could move this down into Particles.create()
      if (max_local_nm < 16*(MAX_PIPELINE+1))
	max_local_nm = 16*(MAX_PIPELINE+1);
#endif
    }
    auto sp = this->create(name, q, m, max_local_np, max_local_nm,
			   sort_interval, sort_out_of_place, vgrid_);
    return this->append(sp);
  }

  static const Convert convert_to_, convert_from_;
  const Convert& convert_to() override { return convert_to_; }
  const Convert& convert_from() override { return convert_from_; }

  const Grid& vgrid() { return *vgrid_; }
  
private:
  Grid* vgrid_;
};

