
#include <psc_balance.h>

#include <psc_fields_as_single.h>
#include <fields.hxx>
#include <bnd.hxx>
#include <inject.hxx>
#include <fields_item.hxx>

#include "cuda_iface.h"

#include <stdlib.h>

void psc_mparticles_cuda_inject(struct psc_mparticles *mprts_base, struct cuda_mparticles_prt *buf,
				uint *buf_n_by_patch); // FIXME

// ======================================================================
// InjectCuda

template<typename Target_t>
struct InjectCuda : InjectBase
{
  using Self = InjectCuda;
  using fields_t = mfields_t::fields_t;
  using Fields = Fields3d<fields_t>;

  InjectCuda(MPI_Comm comm, bool do_inject, int every_step, int tau, int kind_n,
	     Target_t target)
    : InjectBase(do_inject, every_step, tau, kind_n)

  {
    item_n_bnd = psc_bnd_create(comm);
    psc_bnd_set_name(item_n_bnd, "inject_item_n_bnd");
    psc_bnd_set_type(item_n_bnd, "cuda");
    psc_bnd_set_psc(item_n_bnd, ppsc);

    item_n = psc_output_fields_item_create(comm);
    psc_output_fields_item_set_type(item_n, "n_1st_cuda");
    psc_output_fields_item_set_psc_bnd(item_n, item_n_bnd);

    psc_output_fields_item_setup(item_n);
    psc_bnd_setup(item_n_bnd);
  }

  // ----------------------------------------------------------------------
  // dtor

  ~InjectCuda()
  {
    // FIXME, more cleanup needed
  }
  
  // ----------------------------------------------------------------------
  // get_n_in_cell
  //
  // helper function for partition / particle setup FIXME duplicated

  static inline int
  get_n_in_cell(struct psc *psc, struct psc_particle_npt *npt)
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

  static void
  _psc_setup_particle(struct psc *psc, struct cuda_mparticles_prt *cprt,
		      struct psc_particle_npt *npt, int p, double xx[3])
  {
    const Grid_t& grid = psc->grid();
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
  
    assert(npt->kind >= 0 && npt->kind < psc->nr_kinds);
    assert(npt->q == psc->kinds[npt->kind].q);
    assert(npt->m == psc->kinds[npt->kind].m);

    cprt->xi[0] = xx[0] - grid.patches[p].xb[0];
    cprt->xi[1] = xx[1] - grid.patches[p].xb[1];
    cprt->xi[2] = xx[2] - grid.patches[p].xb[2];
    cprt->pxi[0] = pxi;
    cprt->pxi[1] = pyi;
    cprt->pxi[2] = pzi;
    cprt->kind = npt->kind;
    cprt->qni_wni = psc->kinds[npt->kind].q;
  }	      

  // ----------------------------------------------------------------------
  // run

  void run(PscMparticlesBase mprts_base, PscMfieldsBase mflds_base) override
  {
    struct psc *psc = ppsc;

    float fac = 1. / psc->coeff.cori * 
      (every_step * psc->dt / tau) /
      (1. + every_step * psc->dt / tau);

    FieldsItemBase* item = PscFieldsItemBase{item_n}.sub();
    item->run(mflds_base, mprts_base);

    mfields_t mf_n = item->mres().get_as<mfields_t>(kind_n, kind_n+1);

    static struct cuda_mparticles_prt *buf;
    static uint buf_n_alloced;
    if (!buf) {
      buf_n_alloced = 1000;
      buf = (struct cuda_mparticles_prt *) calloc(buf_n_alloced, sizeof(*buf));
    }
    uint buf_n_by_patch[psc->n_patches()];

    uint buf_n = 0;
    psc_foreach_patch(psc, p) {
      buf_n_by_patch[p] = 0;
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

	    if (!target_.is_inside(xx)) {
	      continue;
	    }

	    int n_q_in_cell = 0;
	    for (int kind = 0; kind < nr_pop; kind++) {
	      struct psc_particle_npt npt = {};
	      if (kind < psc->nr_kinds) {
		npt.kind = kind;
		npt.q    = psc->kinds[kind].q;
		npt.m    = psc->kinds[kind].m;
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
		  n_in_cell = get_n_in_cell(psc, &npt);
		}
		n_q_in_cell += npt.q * n_in_cell;
	      } else {
		// FIXME, should handle the case where not the last population is neutralizing
		assert(psc->prm.neutralizing_population == nr_pop - 1);
		n_in_cell = -n_q_in_cell / npt.q;
	      }

	      if (buf_n + n_in_cell > buf_n_alloced) {
		buf_n_alloced = 2 * (buf_n + n_in_cell);
		buf = (struct cuda_mparticles_prt *) realloc(buf, buf_n_alloced * sizeof(*buf));
	      }
	      for (int cnt = 0; cnt < n_in_cell; cnt++) {
		_psc_setup_particle(psc, &buf[buf_n + cnt], &npt, p, xx);
		assert(psc->prm.fractional_n_particles_per_cell);
	      }
	      buf_n += n_in_cell;
	      buf_n_by_patch[p] += n_in_cell;
	    }
	  }
	}
      }
    }

    mf_n.put_as(item->mres(), 0, 0);

    psc_mparticles_cuda_inject(mprts_base.mprts(), buf, buf_n_by_patch);
  }

private:
  Target_t target_;
  struct psc_output_fields_item *item_n;
  struct psc_bnd *item_n_bnd;
  int balance_generation_cnt = {};
};

