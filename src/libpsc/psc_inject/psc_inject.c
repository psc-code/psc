
#include "psc_inject_private.h"

#include <psc_particles_as_single.h> // FIXME
#include <psc_fields_c.h> // FIXME

#include <stdlib.h>

// ======================================================================
// psc_inject

// ----------------------------------------------------------------------
// _psc_inject_create

static void
_psc_inject_create(struct psc_inject *inject)
{
  psc_bnd_set_type(inject->item_n_bnd, "single");
  psc_bnd_set_psc(inject->item_n_bnd, ppsc);

  psc_output_fields_item_set_type(inject->item_n, "n_1st_single");
  psc_output_fields_item_set_psc_bnd(inject->item_n, inject->item_n_bnd);
}

// ----------------------------------------------------------------------
// _psc_inject_setup

static void
_psc_inject_setup(struct psc_inject *inject)
{
  // set up necessary bits for calculating / averaging density moment

  psc_inject_setup_member_objs(inject);

  inject->mflds_n = psc_output_fields_item_create_mfields(inject->item_n);
  psc_mfields_set_name(inject->mflds_n, "mflds_n");
}

// ----------------------------------------------------------------------
// _psc_inject_destroy

static void
_psc_inject_destroy(struct psc_inject *inject)
{
  psc_mfields_destroy(inject->mflds_n);
}

// ----------------------------------------------------------------------
// debug_dump

static void
copy_to_mrc_fld(struct mrc_fld *m3, struct psc_mfields *mflds)
{
  psc_foreach_patch(ppsc, p) {
    struct psc_fields *flds = psc_mfields_get_patch(mflds, p);
    struct mrc_fld_patch *m3p = mrc_fld_patch_get(m3, p);
    mrc_fld_foreach(m3, ix,iy,iz, 0,0) {
      for (int m = 0; m < mflds->nr_fields; m++) {
	MRC_M3(m3p, m, ix,iy,iz) = F3_C(flds, m, ix,iy,iz);
      }
    } mrc_fld_foreach_end;
    mrc_fld_patch_put(m3);
  }
}

static void _mrc_unused
debug_dump(struct mrc_io *io, struct psc_mfields *mflds)
{
  /* if (ppsc->timestep % debug_every_step != 0) { */
  /*   return; */
  /* } */

  struct mrc_fld *mrc_fld = mrc_domain_m3_create(ppsc->mrc_domain);
  mrc_fld_set_name(mrc_fld, psc_mfields_name(mflds));
  mrc_fld_set_param_int(mrc_fld, "nr_ghosts", 2);
  mrc_fld_set_param_int(mrc_fld, "nr_comps", mflds->nr_fields);
  mrc_fld_setup(mrc_fld);
  for (int m = 0; m < mflds->nr_fields; m++) {
    mrc_fld_set_comp_name(mrc_fld, m, psc_mfields_comp_name(mflds, m));
  }
  copy_to_mrc_fld(mrc_fld, mflds);
  mrc_fld_write(mrc_fld, io);
  mrc_fld_destroy(mrc_fld);
}

// ----------------------------------------------------------------------
// calc_n

static void
calc_n(struct psc_inject *inject, struct psc_mparticles *mprts_base,
       struct psc_mfields *mflds_base)
{
  psc_output_fields_item_run(inject->item_n, mflds_base, mprts_base, inject->mflds_n);
#if 0
  static struct mrc_io *io;
  if (!io) {
    io = mrc_io_create(psc_comm(ppsc));
    mrc_io_set_type(io, "xdmf_collective");
    mrc_io_set_param_string(io, "basename", "calc_n");
    mrc_io_set_from_options(io);
    mrc_io_setup(io);
  }

  mrc_io_open(io, "w", psc->timestep, psc->timestep * psc->dt);
  debug_dump(io, inject->mflds_n);
  mrc_io_close(io);
#endif
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
// psc_inject_run
//
// FIXME mostly duplicated from psc_setup_particles

void
psc_inject_run(struct psc_inject *inject, struct psc_mparticles *mprts_base,
	       struct psc_mfields *mflds_base)
{
  struct psc *psc = ppsc;
  
  if (!inject->do_inject ||
      psc->timestep % inject->every_step != 0) {
    return;
  }

  particle_real_t fac = 1. / psc->coeff.cori * 
    (inject->every_step * psc->dt / inject->tau) /
    (1. + inject->every_step * psc->dt / inject->tau);

  calc_n(inject, mprts_base, mflds_base);
  int kind_n = inject->kind_n;
  
  struct psc_mparticles *mprts = psc_mparticles_get_as(mprts_base, PARTICLE_TYPE, 0);
  
  psc_foreach_patch(psc, p) {
    struct psc_particles *prts = psc_mparticles_get_patch(mprts, p);
    struct psc_fields *flds_n = psc_mfields_get_patch(inject->mflds_n, p);
    int *ldims = psc->patch[p].ldims;
    
    int i = prts->n_part;
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
	    psc_ops(psc)->init_npt(psc, kind, xx, &npt);
	    
	    int n_in_cell;
	    if (kind != psc->prm.neutralizing_population) {
	      if (psc->timestep >= 0) {
		npt.n -= F3_C(flds_n, kind_n, jx,jy,jz);
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
	    particles_realloc(prts, i + n_in_cell);
	    for (int cnt = 0; cnt < n_in_cell; cnt++) {
	      particle_t *prt = particles_get_one(prts, i++);
	      
	      _psc_setup_particle(psc, prt, &npt, p, xx);
	      assert(psc->prm.fractional_n_particles_per_cell);
	      prt->qni_wni = psc->kinds[prt->kind].q;
	    }
	  }
	}
      }
    }
    prts->n_part = i;
  }

  psc_mparticles_put_as(mprts, mprts_base, 0);
}

// ----------------------------------------------------------------------
// psc_inject class

struct mrc_class_psc_inject mrc_class_psc_inject = {
  .name             = "psc_inject",
  .size             = sizeof(struct psc_inject),
  .param_descr      = psc_inject_descr,
  .create           = _psc_inject_create,
  .setup            = _psc_inject_setup,
  .destroy          = _psc_inject_destroy,
};

