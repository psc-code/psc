
#ifndef PSC_PARTICLES_VPIC_H
#define PSC_PARTICLES_VPIC_H

#include "psc_particles_private.h"
#include "psc_particles_single.h"

#include "particles.hxx"
#include "particles_traits.hxx"

#include "../libpsc/vpic/vpic_iface.h" // FIXME path

#include "psc_method.h"

struct particle_vpic_t
{
  using real_t = float;
};

struct MparticlesVpic;

// ======================================================================
// ParticlesVpic

struct ParticlesVpic
{
  using real_t = float;
  using Real3 = Vec3<real_t>;
  using Double3 = Vec3<double>;
  
  struct const_accessor
  {
    const_accessor(int n, const ParticlesVpic& prts)
      : n_{n}, prts_{prts}
    {}

    Real3 momentum()   const;// { return prt_.p; }
    real_t w()         const;// { return prt_.w; }
    int kind()         const;// { return prt_.kind; }

    Double3 position() const;
    // {
    //   auto& patch = prts_.grid().patches[prts_.p_];

    //   return patch.xb +
    // 	Double3{prt_.x};
    // }
    
  private:
    int n_;
    const ParticlesVpic& prts_;
  };
  
  ParticlesVpic(MparticlesVpic& mprts)
    : mprts_{mprts}
  {}

  uint size() const;
  const_accessor get(int n) const;
  
private:
  MparticlesVpic& mprts_;
};

// ======================================================================
// MparticlesVpic

struct MparticlesVpic : MparticlesBase
{
  using Self = MparticlesVpic;
  using particle_t = particle_vpic_t; // FIXME, don't have it, but needed here...

  Particles& vmprts_;
  Simulation *sim_;

  // ----------------------------------------------------------------------
  // ctor

  MparticlesVpic(const Grid_t& grid, Grid* vgrid = nullptr)
    : MparticlesBase(grid),
      vmprts_{*new Particles},
      sim_{},
      vgrid_(vgrid)
  {
    assert(grid.n_patches() == 1);
    if (ppsc) {
      psc_method_get_param_ptr(ppsc->method, "sim", (void **) &sim_);
      vgrid_ = sim_->grid_;
    }
  }

  // ----------------------------------------------------------------------
  // get_n_prts

  int get_n_prts() const override
  {
    int n_prts = 0;
    for (auto sp = vmprts_.cbegin(); sp != vmprts_.cend(); ++sp) {
      n_prts += sp->np;
    }
    
    return n_prts;
  }

  // ----------------------------------------------------------------------
  // get_size_all
  
  void get_size_all(uint *n_prts_by_patch) const override
  {
    n_prts_by_patch[0] = get_n_prts();
  }

  // ----------------------------------------------------------------------
  // reserve_all
  //
  // This is a bit iffy, since we don't really want to reallocate stuff here,
  // at least for now, and we wouldn't be able to know how to split this into
  // the different species, anyway.

  void reserve_all(const uint *n_prts_by_patch) override
  {
    for (int p = 0; p < n_patches(); p++) {
      int n_prts = 0, n_prts_alloced = 0;
      for (auto sp = vmprts_.cbegin(); sp != vmprts_.cend(); ++sp) {
	n_prts += sp->np;
	n_prts_alloced += sp->max_np;
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
  
  void resize_all(const uint *n_prts_by_patch) override
  {
    // we can't resize to the numbers given, unless it's "resize to 0", we'll just do nothing
    // The mparticles conversion function should call resize_all() itself first, resizing to
    // 0, and then using push_back, which will increase the count back to the right value
    
    if (n_prts_by_patch[0] == 0) {
      for (auto sp = vmprts_.begin(); sp != vmprts_.end(); ++sp) {
	sp->np = 0;
      }
    } else {
#if 0
      int cur_n_prts_by_patch[n_patches];
      vpic_mparticles_get_size_all(vmprts, n_patches, cur_n_prts_by_patch);
      
      mprintf("vpic_mparticles_resize_all: ignoring %d -> %d\n",
	      cur_n_prts_by_patch[0], n_prts_by_patch[0]);
#endif
    }
  }

  void inject(int p, const particle_inject& prt) override
  {
    particle_inject prt_reweighted = prt;
    //prt_reweighted.w *= dVi;
    inject_reweight(p, prt_reweighted);
  }

  void inject_reweight(int p, const particle_inject& prt) override
  {
    assert(p == 0);
    vmprts_.inject_particle(vmprts_, prt);
  }

  void push_back(const vpic_mparticles_prt *prt)
  {
    for (auto sp = vmprts_.begin(); sp != vmprts_.end(); ++sp) {
      if (sp->id == prt->kind) {
	assert(sp->np < sp->max_np);
	// the below is inject_particle_raw()
	Particles::Particle * RESTRICT p = sp->p + (sp->np++);
	p->dx = prt->dx[0]; p->dy = prt->dx[1]; p->dz = prt->dx[2]; p->i = prt->i;
	p->ux = prt->ux[0]; p->uy = prt->ux[1]; p->uz = prt->ux[2]; p->w = prt->w;
	return;
      }
    }
    mprintf("prt->kind %d not found in species list!\n", prt->kind);
    assert(0);
  }
  
  ParticlesVpic operator[](int p)
  {
    assert(p == 0);
    return ParticlesVpic(*this);
  }

  Particles::Species* define_species(const char *name, double q, double m,
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
    auto sp = vmprts_.create(name, q, m, max_local_np, max_local_nm,
			     sort_interval, sort_out_of_place, vgrid_);
    return vmprts_.append(sp);
  }

  static const Convert convert_to_, convert_from_;
  const Convert& convert_to() override { return convert_to_; }
  const Convert& convert_from() override { return convert_from_; }

private:
  Grid* vgrid_;
};

template<>
struct Mparticles_traits<MparticlesVpic>
{
  static constexpr const char* name = "vpic";
  static MPI_Datatype mpi_dtype() { return MPI_FLOAT; }
};

// ======================================================================

inline uint ParticlesVpic::size() const
{
  return mprts_.get_n_prts();
}

#endif
