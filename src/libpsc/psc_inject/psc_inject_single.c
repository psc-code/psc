
#include <psc_inject_private.h>
#include <psc_balance.h>

#include <psc_particles_as_single.h>
#include <psc_fields_as_c.h>

#include <stdlib.h>

void psc_bnd_check_domain(struct psc_bnd *bnd); // FIXME

// ======================================================================
// psc_inject subclass "single"

// ----------------------------------------------------------------------
// psc_inject_single_create

static void
psc_inject_single_create(struct psc_inject *inject)
{
  // it looks like n_1st_single takes "single" particles, but makes
  // moment fields of type "c", so let's use those "c" fields.
  psc_bnd_set_name(inject->item_n_bnd, "inject_item_n_bnd");
  psc_bnd_set_type(inject->item_n_bnd, "c");
  psc_bnd_set_psc(inject->item_n_bnd, ppsc);

  psc_output_fields_item_set_type(inject->item_n, "n_1st_single");
  psc_output_fields_item_set_psc_bnd(inject->item_n, inject->item_n_bnd);
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
_psc_setup_particle(struct psc *psc, particle_t *prt, struct psc_particle_npt *npt,
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
  
  assert(npt->kind >= 0 && npt->kind < psc->nr_kinds);
  prt->kind = npt->kind;
  assert(npt->q == psc->kinds[prt->kind].q);
  assert(npt->m == psc->kinds[prt->kind].m);
  /* prt->qni = psc->kinds[prt->kind].q; */
  /* prt->mni = psc->kinds[prt->kind].m; */
  prt->xi = xx[0] - psc->patch[p].xb[0];
  prt->yi = xx[1] - psc->patch[p].xb[1];
  prt->zi = xx[2] - psc->patch[p].xb[2];
  prt->pxi = pxi * cos(psc->prm.theta_xz) + pzi * sin(psc->prm.theta_xz);
  prt->pyi = pyi;
  prt->pzi = - pxi * sin(psc->prm.theta_xz) + pzi * cos(psc->prm.theta_xz);
}	      

// ----------------------------------------------------------------------
// psc_inject_single_run

static void
psc_inject_single_run(struct psc_inject *inject, struct psc_mparticles *mprts_base,
		       struct psc_mfields *mflds_base)
{
  struct psc *psc = ppsc;
    
  particle_real_t fac = 1. / psc->coeff.cori * 
    (inject->every_step * psc->dt / inject->tau) /
    (1. + inject->every_step * psc->dt / inject->tau);

  if (psc_balance_generation_cnt != inject->balance_generation_cnt) {
    inject->balance_generation_cnt = psc_balance_generation_cnt;
    psc_bnd_check_domain(inject->item_n_bnd);
  }
  psc_output_fields_item_run(inject->item_n, mflds_base, mprts_base, inject->mflds_n);

  int kind_n = inject->kind_n;
  
  struct psc_mparticles *mprts = psc_mparticles_get_as(mprts_base, PARTICLE_TYPE, 0);
  struct psc_mfields *mflds_n = psc_mfields_get_as(inject->mflds_n, FIELDS_TYPE, kind_n, kind_n+1);
  
  psc_foreach_patch(psc, p) {
    struct psc_fields *flds_n = psc_mfields_get_patch(mflds_n, p);
    int *ldims = psc->patch[p].ldims;
    
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

	  if (!psc_target_is_inside(inject->target, xx)) {
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
	    psc_target_init_npt(inject->target, kind, xx, &npt);
	    
	    int n_in_cell;
	    if (kind != psc->prm.neutralizing_population) {
	      if (psc->timestep >= 0) {
		npt.n -= F3(flds_n, kind_n, jx,jy,jz);
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
	    psc_mparticles_realloc(mprts, p, mparticles_get_n_prts(mprts, p) + n_in_cell);
	    for (int cnt = 0; cnt < n_in_cell; cnt++) {
	      assert(psc->prm.fractional_n_particles_per_cell);
	      particle_t prt;
	      _psc_setup_particle(psc, &prt, &npt, p, xx);
	      prt.qni_wni = psc->kinds[prt.kind].q; // ??? FIXME

	      mparticles_push_back(mprts, p, prt);
	    }
	  }
	}
      }
    }
  }

  psc_mparticles_put_as(mprts, mprts_base, 0);
  psc_mfields_put_as(mflds_n, inject->mflds_n, 0, 0);
}

// ----------------------------------------------------------------------
// psc_inject "single"

struct psc_inject_ops psc_inject_ops_single = {
  .name                = "single",
  .create              = psc_inject_single_create,
  .run                 = psc_inject_single_run,
};

