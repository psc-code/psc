
#pragma once

#include "inject.hxx"
#include "cuda_iface.h"
#include "fields_item_moments_1st_cuda.hxx"

// ======================================================================
// InjectCuda

template<typename _Mparticles, typename _Mfields, typename Target_t, typename dim>
struct InjectCuda : InjectBase
{
  static_assert(_Mparticles::is_cuda::value, "InjectCuda only works with MparticlesCuda");
  
  using Mfields = _Mfields;
  using Mparticles = _Mparticles;
  using fields_t = typename Mfields::fields_t;
  using Fields = Fields3d<fields_t>;
  using real_t = typename Mparticles::real_t;
  using ItemMoment_t = Moment_n_1st_cuda<Mparticles, dim>;
  
  // ----------------------------------------------------------------------
  // ctor
  
  InjectCuda(const Grid_t& grid, int interval, int tau, int kind_n,
	     Target_t target)
    : InjectBase{interval, tau, kind_n},
      target_{target},
      moment_n_{grid, grid.comm()} // FIXME, should just take grid
  {}

  // ----------------------------------------------------------------------
  // dtor

  ~InjectCuda()
  {
    // FIXME, more cleanup needed
  }
  
  // ----------------------------------------------------------------------
  // operator()

  void operator()(Mparticles& mprts)
  {
    const auto& grid = mprts.grid();
    const auto& kinds = grid.kinds;

    SetupParticles<Mparticles> setup_particles;
    // FIXME, this is taken from and kinda specific to psc_flatfoil_yz.cxx,
    // and really shouldn't be replicated so many times anyway
    setup_particles.fractional_n_particles_per_cell = true;
    setup_particles.neutralizing_population = 1;// MY_ELECTRON;
    
    real_t fac = 1. / grid.norm.cori * 
      (interval * grid.dt / tau) / (1. + interval * grid.dt / tau);

    moment_n_.run(mprts);

    MfieldsCuda& mres = moment_n_.result();
    auto& mf_n = mres.get_as<MfieldsSingle>(kind_n, kind_n+1);

    static std::vector<particle_inject> buf2;

    uint buf_n_by_patch[mprts.n_patches()];

    buf2.clear();
    for (int p = 0; p < mprts.n_patches(); p++) {
      buf_n_by_patch[p] = 0;
      Fields N(mf_n[p]);
      const int *ldims = grid.ldims;
    
      for (int jz = 0; jz < ldims[2]; jz++) {
	for (int jy = 0; jy < ldims[1]; jy++) {
	  for (int jx = 0; jx < ldims[0]; jx++) {
	    double xx[3] = {grid.patches[p].x_cc(jx), grid.patches[p].y_cc(jy), grid.patches[p].z_cc(jz)};
	    // FIXME, the issue really is that (2nd order) particle pushers
	    // don't handle the invariant dim right
	    if (grid.isInvar(0) == 1) xx[0] = grid.patches[p].x_nc(jx);
	    if (grid.isInvar(1) == 1) xx[1] = grid.patches[p].y_nc(jy);
	    if (grid.isInvar(2) == 1) xx[2] = grid.patches[p].z_nc(jz);

	    if (!target_.is_inside(xx)) {
	      continue;
	    }

	    int n_q_in_cell = 0;
	    for (int kind = 0; kind < kinds.size(); kind++) {
	      struct psc_particle_npt npt = {};
	      npt.kind = kind;
	      npt.q    = kinds[kind].q;
	      npt.m    = kinds[kind].m;
	      target_.init_npt(kind, xx, &npt);
	    
	      int n_in_cell;
	      if (kind != setup_particles.neutralizing_population) {
		if (grid.timestep() >= 0) {
		  npt.n -= N(kind_n, jx,jy,jz);
		  if (npt.n < 0) {
		    n_in_cell = 0;
		  } else {
		    // this rounds down rather than trying to get fractional particles
		    // statistically right...
		    n_in_cell = npt.n *fac;		}
		} else {
		  n_in_cell = setup_particles.get_n_in_cell(grid, &npt);
		}
		n_q_in_cell += npt.q * n_in_cell;
	      } else {
		// FIXME, should handle the case where not the last population is neutralizing
		assert(setup_particles.neutralizing_population == kinds.size() - 1);
		n_in_cell = -n_q_in_cell / npt.q;
	      }

	      for (int cnt = 0; cnt < n_in_cell; cnt++) {
		assert(setup_particles.fractional_n_particles_per_cell);
		real_t wni = 1.; // ??? FIXME
		auto prt = particle_inject{{}, {}, wni, npt.kind};
		setup_particles.setup_particle(grid, &prt, &npt, p, xx);
		buf2.push_back(prt);
	      }
	      buf_n_by_patch[p] += n_in_cell;
	    }
	  }
	}
      }
    }

    mres.put_as(mf_n, 0, 0);

    auto cur = buf2.cbegin();
    for (int p = 0; p < mprts.n_patches(); p++) {
      auto injector = mprts[p].injector();
      for (int n = 0; n < buf_n_by_patch[p]; n++) {
	injector(*cur++);
      }
    }
    mprintf("**** Inject: %ld particles added\n", buf2.size());
  }

  // ----------------------------------------------------------------------
  // run
  
  void run(MparticlesBase& mprts_base, MfieldsBase& mflds_base) override
  {
    auto& mprts = mprts_base.get_as<Mparticles>();
    (*this)(mprts);
    mprts_base.put_as(mprts);
  }
  
private:
  Target_t target_;
  ItemMoment_t moment_n_;
};

