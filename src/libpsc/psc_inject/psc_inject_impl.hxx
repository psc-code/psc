
#include <psc_balance.h>

#include <inject.hxx>
#include <fields.hxx>
#include <bnd.hxx>
#include <fields_item.hxx>
#include "../libpsc/psc_output_fields/fields_item_moments_1st.hxx"

#include <stdlib.h>
#include <string>

// ======================================================================
// Inject_

template<typename MP, typename MF, typename Target_t>
struct Inject_ : InjectBase
{
  using Self = Inject_<MP, MF, Target_t>;
  using Mfields = MF;
  using Mparticles = MP;
  using fields_t = typename Mfields::fields_t;
  using Fields = Fields3d<fields_t>;
  using real_t = typename Mparticles::real_t;
  using particle_t = typename Mparticles::particle_t;
  using mparticles_t = PscMparticles<Mparticles>;
  using mfields_t = PscMfields<Mfields>;
  using ItemMoment_t = ItemMomentLoopPatches<Moment_n_1st<mparticles_t, mfields_t>>;
  
  // ----------------------------------------------------------------------
  // ctor
  
  Inject_(MPI_Comm comm, bool do_inject, int every_step, int tau, int kind_n,
	  Target_t target)
    : InjectBase{do_inject, every_step, tau, kind_n},
      target_{target}
  {
    // it looks like n_1st_sub takes "sub" particles, but makes
    // moment fields of type "c", so let's use those "c" fields.
    item_n_bnd = psc_bnd_create(comm);
    psc_bnd_set_name(item_n_bnd, "inject_item_n_bnd");
    psc_bnd_set_type(item_n_bnd, "c");
    psc_bnd_set_psc(item_n_bnd, ppsc);

    psc_bnd_setup(item_n_bnd);
    moment_n_.reset(new ItemMoment_t{comm, PscBndBase{item_n_bnd}});
  }

  // ----------------------------------------------------------------------
  // dtor

  ~Inject_()
  {
    // FIXME, more cleanup needed
  }
  
  // ----------------------------------------------------------------------
  // run

  void run(PscMparticlesBase mprts_base, PscMfieldsBase mflds_base) override
  {
    mparticles_t mprts = mprts_base.get_as<mparticles_t>();
    (*this)(*mprts.sub());
    mprts.put_as(mprts_base);
  }
  
  void operator()(Mparticles& mprts)
  {
    struct psc *psc = ppsc;
    const auto& grid = mprts.grid();
    const auto& kinds = grid.kinds;
    
    real_t fac = 1. / psc->coeff.cori * 
      (every_step * psc->dt / tau) /
      (1. + every_step * psc->dt / tau);

    moment_n_->run(mprts);
    auto mf_n = moment_n_->mres();

    psc_foreach_patch(psc, p) {
      Fields N(mf_n[p]);
      const int *ldims = psc->grid().ldims;
    
      int nr_pop = psc->prm.nr_populations;
      for (int jz = 0; jz < ldims[2]; jz++) {
	for (int jy = 0; jy < ldims[1]; jy++) {
	  for (int jx = 0; jx < ldims[0]; jx++) {
	    double xx[3] = { .5 * (CRDX(p, jx) + CRDX(p, jx+1)),
			     .5 * (CRDY(p, jy) + CRDY(p, jy+1)),
			     .5 * (CRDZ(p, jz) + CRDZ(p, jz+1)) };
	    // FIXME, the issue really is that (2nd order) particle pushers
	    // don't handle the invariant dim right
	    if (grid.gdims[0] == 1) xx[0] = CRDX(p, jx);
	    if (grid.gdims[1] == 1) xx[1] = CRDY(p, jy);
	    if (grid.gdims[2] == 1) xx[2] = CRDZ(p, jz);

	    if (!target_.is_inside(xx)) {
	      continue;
	    }

	    int n_q_in_cell = 0;
	    for (int kind = 0; kind < nr_pop; kind++) {
	      struct psc_particle_npt npt = {};
	      if (kind < kinds.size()) {
		npt.kind = kind;
		npt.q    = kinds[kind].q;
		npt.m    = kinds[kind].m;
	      };
	      target_.init_npt(kind, xx, &npt);
	    
	      int n_in_cell;
	      if (kind != psc->prm.neutralizing_population) {
		if (psc->timestep >= 0) {
		  npt.n -= N(kind_n, jx,jy,jz);
		  if (npt.n < 0) {
		    n_in_cell = 0;
		  } else {
		    // this rounds down rather than trying to get fractional particles
		    // statistically right...
		    n_in_cell = npt.n *fac;		}
		} else {
		  n_in_cell = SetupParticles<Mparticles>::get_n_in_cell(psc, &npt);
		}
		n_q_in_cell += npt.q * n_in_cell;
	      } else {
		// FIXME, should handle the case where not the last population is neutralizing
		assert(psc->prm.neutralizing_population == nr_pop - 1);
		n_in_cell = -n_q_in_cell / npt.q;
	      }
	      for (int cnt = 0; cnt < n_in_cell; cnt++) {
		assert(psc->prm.fractional_n_particles_per_cell);
		particle_t prt;
		SetupParticles<Mparticles>::setup_particle(psc, &prt, &npt, p, xx);
		prt.qni_wni_ = kinds[prt.kind_].q; // ??? FIXME

		mprts[p].push_back(prt);
	      }
	    }
	  }
	}
      }
    }
  }

private:
  Target_t target_;
  std::unique_ptr<ItemMoment_t> moment_n_;
  struct psc_bnd *item_n_bnd;
  int balance_generation_cnt = {};
};

