
#include <psc.h>
#include <psc_push_particles.h>
#include <psc_moments.h>
#include <psc_push_fields.h>
#include <psc_sort.h>
#include <psc_balance.h>
#include <psc_particles_as_c.h>

#include <mrc_params.h>
#include <mrc_profile.h>

#include <math.h>
#include <time.h>

// ======================================================================

#include "psc_output_fields_item_private.h"

#include "psc_bnd.h"
#include "psc_fields_as_c.h"

static void
do_1st_calc_kh(int p, fields_t *pf, particles_t *pp)
{
  fields_c_real_t fnqs = sqr(ppsc->coeff.alpha) * ppsc->coeff.cori / ppsc->coeff.eta;
  fields_c_real_t dxi = 1.f / ppsc->dx[0];
  fields_c_real_t dyi = 1.f / ppsc->dx[1];
  fields_c_real_t dzi = 1.f / ppsc->dx[2];

  struct psc_patch *patch = &ppsc->patch[p];
  for (int n = 0; n < pp->n_part; n++) {
    particle_t *part = particles_get_one(pp, n);
      
    fields_c_real_t u = (part->xi - patch->xb[0]) * dxi - .5;
    fields_c_real_t v = (part->yi - patch->xb[1]) * dyi - .5;
    fields_c_real_t w = (part->zi - patch->xb[2]) * dzi - .5;
    int j1 = particle_real_fint(u);
    int j2 = particle_real_fint(v);
    int j3 = particle_real_fint(w);
    fields_c_real_t h1 = u-j1;
    fields_c_real_t h2 = v-j2;
    fields_c_real_t h3 = w-j3;
      
    fields_c_real_t g0x=1.f - h1;
    fields_c_real_t g0y=1.f - h2;
    fields_c_real_t g0z=1.f - h3;
    fields_c_real_t g1x=h1;
    fields_c_real_t g1y=h2;
    fields_c_real_t g1z=h3;
      
    if (ppsc->domain.gdims[0] == 1) {
      j1 = 0; g0x = 1.; g1x = 0.;
    }
    if (ppsc->domain.gdims[1] == 1) {
      j2 = 0; g0y = 1.; g1y = 0.;
    }
    if (ppsc->domain.gdims[2] == 1) {
      j3 = 0; g0z = 1.; g1z = 0.;
    }

    assert(j1 >= -1 && j1 < patch->ldims[0]);
    assert(j2 >= -1 && j2 < patch->ldims[1]);
    assert(j3 >= -1 && j3 < patch->ldims[2]);
      
    fields_c_real_t fnq;
    int m;
    if (part->qni < 0.) {
      fnq = part->qni * part->wni * fnqs;
      if (part->wni > 1.) {
	m = 1;
      } else {
	m = 0;
      }
    } else if (part->qni > 0.) {
      fnq = part->qni * part->wni * fnqs;
      if (part->wni > 1.) {
	m = 3;
      } else {
	m = 2;
      }
    } else {
      assert(0);
    }
    F3(pf, m, j1  ,j2  ,j3  ) += fnq*g0x*g0y*g0z;
    F3(pf, m, j1+1,j2  ,j3  ) += fnq*g1x*g0y*g0z;
    F3(pf, m, j1  ,j2+1,j3  ) += fnq*g0x*g1y*g0z;
    F3(pf, m, j1+1,j2+1,j3  ) += fnq*g1x*g1y*g0z;
    F3(pf, m, j1  ,j2  ,j3+1) += fnq*g0x*g0y*g1z;
    F3(pf, m, j1+1,j2  ,j3+1) += fnq*g1x*g0y*g1z;
    F3(pf, m, j1  ,j2+1,j3+1) += fnq*g0x*g1y*g1z;
    F3(pf, m, j1+1,j2+1,j3+1) += fnq*g1x*g1y*g1z;
  }
}

static void
calc_kh(struct psc_output_fields_item *item, mfields_base_t *flds,
	mparticles_base_t *particles_base, mfields_c_t *res)
{
  static int pr;
  if (!pr) {
    pr = prof_register("c_moments_kh", 1., 0, 0);
  }

  mparticles_t *particles = psc_mparticles_get_cf(particles_base, 0);

  prof_start(pr);
  for (int m = 0; m < 4; m++) {
    psc_mfields_zero(res, m);
  }
  
  psc_foreach_patch(ppsc, p) {
    do_1st_calc_kh(p, psc_mfields_get_patch_c(res, p),
		   psc_mparticles_get_patch(particles, p));
  }
  prof_stop(pr);

  psc_mparticles_put_cf(particles, particles_base); // FIXME, don't need copy-back

  psc_bnd_add_ghosts(item->bnd, res, 0, 4);
  // FIXME add_ghosts_boundary(res, 0, 4);
}

static struct psc_output_fields_item_ops psc_output_fields_item_kh_ops = {
  .name      = "kh",
  .nr_comp   = 4,
  .fld_names = { "ne1", "ne2", "ni1", "ni2" },
  .run       = calc_kh,
};

// ======================================================================

struct psc_kh {
  // parameters
  double beta;
  double theta_B, theta_V;
  double delta;
  double mi_over_me;
  double wpe_over_wce;
  double Ti_over_Te;

  // calculated from the above
  double B0;
  double v0z;
  double Te, Ti;
};

#define to_psc_kh(psc) mrc_to_subobj(psc, struct psc_kh)

#define VAR(x) (void *)offsetof(struct psc_kh, x)
static struct param psc_kh_descr[] = {
  { "theta_B"       , VAR(theta_B)         , PARAM_DOUBLE(M_PI/2. - .05) },
  { "theta_V"       , VAR(theta_V)         , PARAM_DOUBLE(M_PI/2. - .05) },
  { "delta"         , VAR(delta)           , PARAM_DOUBLE(2.)            },
  { "beta"          , VAR(beta)            , PARAM_DOUBLE(.5)            },
  { "mi_over_me"    , VAR(mi_over_me)      , PARAM_DOUBLE(5.)            },
  { "wpe_over_wce"  , VAR(wpe_over_wce)    , PARAM_DOUBLE(2.)            },
  { "Ti_over_Te"    , VAR(Ti_over_Te)      , PARAM_DOUBLE(1.)            },
  {},
};
#undef VAR

// ----------------------------------------------------------------------
// psc_kh_create

static void
psc_kh_create(struct psc *psc)
{
  // new defaults (dimensionless) for this case
  psc->prm.qq = 1.;
  psc->prm.mm = 1.;
  psc->prm.tt = 1.;
  psc->prm.cc = 1.;
  psc->prm.eps0 = 1.;

  psc->prm.nmax = 16000;
  psc->prm.cpum = 5*24.0*60*60;
  psc->prm.lw = 2.*M_PI;
  psc->prm.i0 = 0.;
  psc->prm.n0 = 1.;
  psc->prm.e0 = 1.;

  psc->prm.nicell = 50;
  psc->prm.gdims_in_terms_of_cells = true;
  psc->prm.cfl = 0.98;

  psc->domain.length[0] = 1.; // no x-dependence
  psc->domain.length[1] = 30.;
  psc->domain.length[2] = 30.;

  psc->domain.gdims[0] = 1;
  psc->domain.gdims[1] = 160;
  psc->domain.gdims[2] = 160;

  psc->domain.bnd_fld_lo[0] = BND_FLD_PERIODIC;
  psc->domain.bnd_fld_hi[0] = BND_FLD_PERIODIC;
  psc->domain.bnd_fld_lo[1] = BND_FLD_CONDUCTING_WALL;
  psc->domain.bnd_fld_hi[1] = BND_FLD_CONDUCTING_WALL;
  psc->domain.bnd_fld_lo[2] = BND_FLD_PERIODIC;
  psc->domain.bnd_fld_hi[2] = BND_FLD_PERIODIC;
  psc->domain.bnd_part[0] = BND_PART_PERIODIC;
  psc->domain.bnd_part[1] = BND_PART_REFLECTING;
  psc->domain.bnd_part[2] = BND_PART_PERIODIC;

  // FIXME: can only use 1st order pushers with current conducting wall b.c.
  psc_push_particles_set_type(psc->push_particles, "1vb");
  psc_moments_set_type(psc->moments, "1st_cc");
}

// ----------------------------------------------------------------------
// psc_kh_setup
//
// the parameters are now set, calculate quantities to initialize fields,
// particles

static void
psc_kh_setup(struct psc *psc)
{
  struct psc_kh *kh = to_psc_kh(psc);

  double me = 1.;
  double B0 = sqrt(me) / (kh->wpe_over_wce);
  double vAe = B0 / sqrt(me);
  double vAe_plane = vAe * cos(kh->theta_V);
  double v0 = 2 * vAe_plane;
  double v0z = v0 * cos(kh->theta_V - kh->theta_B);
  double Te = kh->beta * (1. / (1. + kh->Ti_over_Te)) * sqr(B0) / 2.;
  double Ti = kh->beta * (1. / (1. + 1./kh->Ti_over_Te)) * sqr(B0) / 2.;

  kh->B0 = B0;
  kh->v0z = v0z;
  kh->Te = Te;
  kh->Ti = Ti;

  psc_setup_default(psc);
}

// ----------------------------------------------------------------------
// psc_kh_init_field

static double
psc_kh_init_field(struct psc *psc, double x[3], int m)
{
  struct psc_kh *kh = to_psc_kh(psc);

  double yl = psc->domain.length[1];
  double vz = kh->v0z * tanh((x[1] - .5 * yl) / kh->delta);

  switch (m) {
  case HX: return kh->B0 * sin(kh->theta_B);
  case HZ: return kh->B0 * cos(kh->theta_B);
  case EY: return -vz * kh->B0 * sin(kh->theta_B);
  default: return 0.;
  }
}

// ----------------------------------------------------------------------
// psc_kh_init_npt

static void
psc_kh_init_npt(struct psc *psc, int kind, double x[3],
		struct psc_particle_npt *npt, double *wni)
{
  struct psc_kh *kh = to_psc_kh(psc);

  double yl = psc->domain.length[1];
  double vz = kh->v0z * tanh((x[1] - .5 * yl) / kh->delta);

  npt->n = 1.;
  npt->p[2] = vz;

  if (vz < 0) {
    *wni = 1.;
  } else {
    *wni = 1. + 1e-6;
  }

  switch (kind) {
  case 0: // electrons
    npt->q = -1.;
    npt->m = 1.;
    npt->T[0] = kh->Te;
    npt->T[1] = kh->Te;
    npt->T[2] = kh->Te;
    break;
  case 1: // ions
    npt->q = 1.;
    npt->m = kh->mi_over_me;
    npt->T[0] = kh->Ti;
    npt->T[1] = kh->Ti;
    npt->T[2] = kh->Ti;
    break;
  default:
    assert(0);
  }
}

// ----------------------------------------------------------------------
// psc_kh_setup_particles

void
psc_kh_setup_particles(struct psc *psc, int *nr_particles_by_patch,
		       bool count_only)
{
  double beta = psc->coeff.beta;

  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  if (psc->prm.seed_by_time) {
    srandom(10*rank + time(NULL));
  } else {
    srandom(rank);
  }

  mparticles_t *particles = NULL;
  if (!count_only) {
    particles = psc_mparticles_get_cf(psc->particles, MP_DONT_COPY);
  }

  psc_foreach_patch(psc, p) {
    particles_t *pp = NULL;
    if (!count_only) {
      pp = psc_mparticles_get_patch(particles, p);
    }

    int *ldims = psc->patch[p].ldims;
    int i = 0;
    for (int kind = 0; kind < psc->prm.nr_kinds; kind++) {
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

	    struct psc_particle_npt npt = {}; // init to all zero
	    double wni;
	    psc_kh_init_npt(psc, kind, xx, &npt, &wni);
	    
	    int n_in_cell = npt.n * psc->prm.nicell + .5;
	    if (count_only) {
	      i += n_in_cell;
	    }
	    if (!count_only) {
	      for (int cnt = 0; cnt < n_in_cell; cnt++) {
		particle_t *p = particles_get_one(pp, i++);
		
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
		
		float px =
		  sqrtf(-2.f*npt.T[0]/npt.m*sqr(beta)*logf(1.0-ran1)) * cosf(2.f*M_PI*ran2)
		  + npt.p[0];
		float py =
		  sqrtf(-2.f*npt.T[1]/npt.m*sqr(beta)*logf(1.0-ran3)) * cosf(2.f*M_PI*ran4)
		  + npt.p[1];
		float pz =
		  sqrtf(-2.f*npt.T[2]/npt.m*sqr(beta)*logf(1.0-ran5)) * cosf(2.f*M_PI*ran6)
		  + npt.p[2];
		
		p->xi = xx[0];
		p->yi = xx[1];
		p->zi = xx[2];
		p->pxi = px;
		p->pyi = py;
		p->pzi = pz;
		p->qni = npt.q;
		p->mni = npt.m;
		p->wni = wni;
	      }
	    }
	  }
	}
      }
    }
    if (count_only) {
      nr_particles_by_patch[p] = i;
    } else {
      pp->n_part = i;
      assert(pp->n_part == nr_particles_by_patch[p]);
    }
  }
  if (!count_only) {
    psc_mparticles_put_cf(particles, psc->particles);
  }
}

// ======================================================================
// psc_kh_ops

struct psc_ops psc_kh_ops = {
  .name             = "kh",
  .size             = sizeof(struct psc_kh),
  .param_descr      = psc_kh_descr,
  .create           = psc_kh_create,
  .setup            = psc_kh_setup,
  .init_field       = psc_kh_init_field,
  .setup_particles  = psc_kh_setup_particles,
};

// ======================================================================
// main

int
main(int argc, char **argv)
{
  mrc_class_register_subclass(&mrc_class_psc_output_fields_item, &psc_output_fields_item_kh_ops);

  return psc_main(&argc, &argv, &psc_kh_ops);
}
