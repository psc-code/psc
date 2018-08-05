
#ifndef PSC_PARTICLES_VPIC_H
#define PSC_PARTICLES_VPIC_H

#include "psc_particles_private.h"
#include "psc_particles_single.h"

#include "particles.hxx"
#include "particles_traits.hxx"

#include "../libpsc/vpic/vpic_iface.h" // FIXME path

#include "psc_method.h"

void vpic_mparticles_get_size_all(Particles *vmprts, int n_patches,
				  uint *n_prts_by_patch);

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

  Particles *vmprts;
  Simulation *sim;

  MparticlesVpic(const Grid_t& grid)
    : MparticlesBase(grid)
  {
    vmprts = new Particles;
    psc_method_get_param_ptr(ppsc->method, "sim", (void **) &sim);
  }

  static mrc_obj_method methods[];

  int get_n_prts() const override
  {
    return sim->mprts_get_nr_particles(*vmprts);
  }
  
  void get_size_all(uint *n_prts_by_patch) const override
  {
    vpic_mparticles_get_size_all(vmprts, 1, n_prts_by_patch);
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
      for (auto sp = vmprts->cbegin(); sp != vmprts->cend(); ++sp) {
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
      for (auto sp = vmprts->begin(); sp != vmprts->end(); ++sp) {
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

  void inject_reweight(int p, const particle_inject& prt) override
  {
    assert(p == 0);
    static_cast<ParticlesOps*>(sim)->inject_particle(*vmprts, *sim->accumulator_, *sim->field_array_, &prt);
  }

  void push_back(const struct vpic_mparticles_prt *prt);

  ParticlesVpic operator[](int p)
  {
    assert(p == 0);
    return ParticlesVpic(*this);
  }

  static const Convert convert_to_, convert_from_;
  const Convert& convert_to() override { return convert_to_; }
  const Convert& convert_from() override { return convert_from_; }
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
