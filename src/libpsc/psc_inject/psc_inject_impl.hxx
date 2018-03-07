
#include <psc_inject_private.h>
#include <psc_balance.h>

#include <inject.hxx>
#include <fields.hxx>
#include <bnd.hxx>

#include <stdlib.h>
#include <string>

// ======================================================================
// Inject_

template<typename MP, typename MF>
struct Inject_ : InjectBase
{
  using Self = Inject_<MP, MF>;
  using Mfields = MF;
  using Mparticles = MP;
  using fields_t = typename Mfields::fields_t;
  using Fields = Fields3d<fields_t>;
  using real_t = typename Mparticles::real_t;
  using particle_t = typename Mparticles::particle_t;
  using mparticles_t = PscMparticles<Mparticles>;
  using mfields_t = PscMfields<Mfields>;
  
  // ----------------------------------------------------------------------
  // ctor
  
  Inject_(MPI_Comm comm, bool do_inject, int every_step, int tau, int kind_n,
	  psc_target* target)
    : InjectBase(do_inject, every_step, tau, kind_n, target)
  {
    // it looks like n_1st_sub takes "sub" particles, but makes
    // moment fields of type "c", so let's use those "c" fields.
    item_n_bnd = psc_bnd_create(comm);
    psc_bnd_set_name(item_n_bnd, "inject_item_n_bnd");
    psc_bnd_set_type(item_n_bnd, "c");
    psc_bnd_set_psc(item_n_bnd, ppsc);

    item_n = psc_output_fields_item_create(comm);
    auto name = std::string("n_1st_") + mparticles_traits<mparticles_t>::name;
    psc_output_fields_item_set_type(item_n, name.c_str());
    psc_output_fields_item_set_psc_bnd(item_n, item_n_bnd);

    mflds_n = psc_output_fields_item_create_mfields(item_n);
    psc_output_fields_item_setup(item_n);
    psc_bnd_setup(item_n_bnd);

    psc_mfields_set_name(mflds_n, "mflds_n");
    psc_mfields_list_add(&psc_mfields_base_list, &mflds_n);
  }

  // ----------------------------------------------------------------------
  // dtor

  ~Inject_()
  {
    psc_mfields_destroy(mflds_n);
    // FIXME, more cleanup needed
  }
  
  // ----------------------------------------------------------------------
  // get_n_in_cell
  //
  // helper function for partition / particle setup FIXME duplicated

  static int get_n_in_cell(struct psc *psc, struct psc_particle_npt *npt)
  {
    if (psc->prm.const_num_particles_per_cell) {
      return psc->prm.nicell;
    }
    if (npt->particles_per_cell) {
      return npt->n * npt->particles_per_cell + .5;
    }
    if (psc->prm.fractional_n_particles_per_cell) {
      int n_prts = npt->n / psc->coeff.cori;
      float rmndr = npt->n / psc->coeff.cori - n_prts;
      float ran = random() / ((float) RAND_MAX + 1);
      if (ran < rmndr) {
	n_prts++;
      }
      return n_prts;
    }
    return npt->n / psc->coeff.cori + .5;
  }

  // FIXME duplicated

  static void _psc_setup_particle(struct psc *psc, particle_t *prt, struct psc_particle_npt *npt,
				  int p, double xx[3])
  {
    double beta = psc->coeff.beta;

    float ran1, ran2, ran3, ran4, ran5, ran6;
    do {
      ran1 = random() / ((float) RAND_MAX + 1);
      ran2 = random() / ((float) RAND_MAX + 1);
      ran3 = random() / ((float) RAND_MAX + 1);
      ran4 = random() / ((float) RAND_MAX + 1);
      ran5 = random() / ((float) RAND_MAX + 1);
      ran6 = random() / ((float) RAND_MAX + 1);
    } while (ran1 >= 1.f || ran2 >= 1.f || ran3 >= 1.f ||
	     ran4 >= 1.f || ran5 >= 1.f || ran6 >= 1.f);
	      
    double pxi = npt->p[0] +
      sqrtf(-2.f*npt->T[0]/npt->m*sqr(beta)*logf(1.0-ran1)) * cosf(2.f*M_PI*ran2);
    double pyi = npt->p[1] +
      sqrtf(-2.f*npt->T[1]/npt->m*sqr(beta)*logf(1.0-ran3)) * cosf(2.f*M_PI*ran4);
    double pzi = npt->p[2] +
      sqrtf(-2.f*npt->T[2]/npt->m*sqr(beta)*logf(1.0-ran5)) * cosf(2.f*M_PI*ran6);

    if (psc->prm.initial_momentum_gamma_correction) {
      double gam;
      if (sqr(pxi) + sqr(pyi) + sqr(pzi) < 1.) {
	gam = 1. / sqrt(1. - sqr(pxi) - sqr(pyi) - sqr(pzi));
	pxi *= gam;
	pyi *= gam;
	pzi *= gam;
      }
    }

    const Grid_t& grid = psc->grid();
    assert(npt->kind >= 0 && npt->kind < psc->nr_kinds);
    prt->kind_ = npt->kind;
    assert(npt->q == psc->kinds[prt->kind_].q);
    assert(npt->m == psc->kinds[prt->kind_].m);
    /* prt->qni = psc->kinds[prt->kind].q; */
    /* prt->mni = psc->kinds[prt->kind].m; */
    prt->xi = xx[0] - grid.patches[p].xb[0];
    prt->yi = xx[1] - grid.patches[p].xb[1];
    prt->zi = xx[2] - grid.patches[p].xb[2];
    prt->pxi = pxi * cos(psc->prm.theta_xz) + pzi * sin(psc->prm.theta_xz);
    prt->pyi = pyi;
    prt->pzi = - pxi * sin(psc->prm.theta_xz) + pzi * cos(psc->prm.theta_xz);
  }	      

  // ----------------------------------------------------------------------
  // run

  void run(PscMparticlesBase mprts_base, PscMfieldsBase mflds_base) override
  {
    struct psc *psc = ppsc;
    
    real_t fac = 1. / psc->coeff.cori * 
      (every_step * psc->dt / tau) /
      (1. + every_step * psc->dt / tau);

    if (psc_balance_generation_cnt != balance_generation_cnt) {
      balance_generation_cnt = psc_balance_generation_cnt;
      auto bnd = PscBndBase(item_n_bnd);
      bnd.reset();
    }
    psc_output_fields_item_run(item_n, mflds_base.mflds(), mprts_base.mprts(), mflds_n);

    mparticles_t mprts = mprts_base.get_as<mparticles_t>();
    mfields_t mf_n = mflds_n->get_as<mfields_t>(kind_n, kind_n+1);

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
	    if (psc->domain.gdims[0] == 1) xx[0] = CRDX(p, jx);
	    if (psc->domain.gdims[1] == 1) xx[1] = CRDY(p, jy);
	    if (psc->domain.gdims[2] == 1) xx[2] = CRDZ(p, jz);

	    if (!psc_target_is_inside(target, xx)) {
	      continue;
	    }

	    int n_q_in_cell = 0;
	    for (int kind = 0; kind < nr_pop; kind++) {
	      struct psc_particle_npt npt = {};
	      if (kind < psc->nr_kinds) {
		npt.kind = kind;
		npt.q    = psc->kinds[kind].q;
		npt.m    = psc->kinds[kind].m;
		npt.n    = psc->kinds[kind].n;
		npt.T[0] = psc->kinds[kind].T;
		npt.T[1] = psc->kinds[kind].T;
		npt.T[2] = psc->kinds[kind].T;
	      };
	      psc_target_init_npt(target, kind, xx, &npt);
	    
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
		  n_in_cell = get_n_in_cell(psc, &npt);
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
		_psc_setup_particle(psc, &prt, &npt, p, xx);
		prt.qni_wni_ = psc->kinds[prt.kind_].q; // ??? FIXME

		mprts[p].push_back(prt);
	      }
	    }
	  }
	}
      }
    }

    mprts.put_as(mprts_base);
    mf_n.put_as(mflds_n, 0, 0);
  }

};

