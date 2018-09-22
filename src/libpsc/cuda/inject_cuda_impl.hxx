
#pragma once

#include "inject.hxx"
#include "cuda_iface.h"
#include "fields_item_moments_1st_cuda.hxx"

// ======================================================================
// InjectCuda

template<typename BS, typename dim, typename Target_t>
struct InjectCuda : InjectBase
{
  using Self = InjectCuda;
  using fields_t = MfieldsSingle::fields_t;
  using Fields = Fields3d<fields_t>;
  using ItemMoment_t = Moment_n_1st_cuda<BS, dim>;
  
  InjectCuda(const Grid_t& grid, int interval, int tau, int kind_n, Target_t target)
    : InjectBase(interval, tau, kind_n),
      target_{target},
      moment_n_{grid, grid.comm()}, // FIXME, should just take grid
      grid_{grid}
  {}

  // ----------------------------------------------------------------------
  // dtor

  ~InjectCuda()
  {
    // FIXME, more cleanup needed
  }
  
  // ----------------------------------------------------------------------
  // operator()

  void operator()(MparticlesCuda<BS>& mprts)
  {
    const auto& grid = mprts.grid();
    const auto& kinds = grid_.kinds;

    SetupParticles<MparticlesCuda<BS>> setup_particles;

    float fac = 1. / grid_.norm.cori * 
      (interval * grid_.dt / tau) / (1. + interval * grid_.dt / tau);

    moment_n_.run(mprts);

    MfieldsCuda& mres = moment_n_.result();
    auto& mf_n = mres.get_as<MfieldsSingle>(kind_n, kind_n+1);

    static std::vector<cuda_mparticles_prt> buf;

    uint buf_n_by_patch[grid_.n_patches()];

    buf.clear();
    for (int p = 0; p < grid_.n_patches(); p++) {
      buf_n_by_patch[p] = 0;
      Fields N(mf_n[p]);
      const int *ldims = grid_.ldims;
    
      for (int jz = 0; jz < ldims[2]; jz++) {
	for (int jy = 0; jy < ldims[1]; jy++) {
	  for (int jx = 0; jx < ldims[0]; jx++) {
	    double xx[3] = {grid_.patches[p].x_cc(jx), grid_.patches[p].y_cc(jy), grid_.patches[p].z_cc(jz)};
	    // FIXME, the issue really is that (2nd order) particle pushers
	    // don't handle the invariant dim right
	    if (grid_.isInvar(0) == 1) xx[0] = grid_.patches[p].x_nc(jx);
	    if (grid_.isInvar(1) == 1) xx[1] = grid_.patches[p].y_nc(jy);
	    if (grid_.isInvar(2) == 1) xx[2] = grid_.patches[p].z_nc(jz);

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
	      if (kind != neutralizing_population) {
		if (grid_.timestep() >= 0) {
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
		assert(neutralizing_population == kinds.size() - 1);
		n_in_cell = -n_q_in_cell / npt.q;
	      }

	      for (int cnt = 0; cnt < n_in_cell; cnt++) {
		cuda_mparticles_prt prt;
		setup_particles.setup_particle(grid, &prt, &npt, p, xx);
		buf.push_back(prt);
		assert(fractional_n_particles_per_cell);
	      }
	      buf_n_by_patch[p] += n_in_cell;
	    }
	  }
	}
      }
    }

    mres.put_as(mf_n, 0, 0);

    mprts.inject_buf(buf.data(), buf_n_by_patch);
    mprintf("**** Inject: %ld particles added\n", buf.size());
  }

  // ----------------------------------------------------------------------
  // run
  
  void run(MparticlesBase& mprts_base, MfieldsBase& mflds_base) override
  {
    auto& mprts = mprts_base.get_as<MparticlesCuda<BS>>();
    (*this)(mprts);
    mprts_base.put_as(mprts);
  }
  
private:
  Target_t target_;
  ItemMoment_t moment_n_;
  const Grid_t& grid_;

  // FIXME
  bool fractional_n_particles_per_cell = { true };
  bool initial_momentum_gamma_correction = { false };
  int neutralizing_population = { 1 };
};

