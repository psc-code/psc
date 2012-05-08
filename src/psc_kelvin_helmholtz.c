
#include <psc.h>
#include <psc_push_particles.h>
#include <psc_moments.h>
#include <psc_push_fields.h>
#include <psc_sort.h>
#include <psc_balance.h>
#include <psc_particles_as_single.h>

#include <mrc_params.h>
#include <mrc_profile.h>

#include <math.h>
#include <time.h>

// ======================================================================

enum {
  KH_ELECTRON1,
  KH_ELECTRON2,
  KH_ION1,
  KH_ION2,
  NR_KH_KINDS,
};

int
particle_single_kind(particle_single_t *prt)
{
  if (particle_qni(prt) < 0.) {
    if (particle_wni(prt) == 1.) {
      return KH_ELECTRON1;
    } else {
      return KH_ELECTRON2;
    }
  } else if (particle_qni(prt) > 0.) {
    if (particle_wni(prt) == 1.) {
      return KH_ION1;
    } else {
      return KH_ION2;
    }
  } else {
    assert(0);
  }
}

// ======================================================================

#include "psc_output_fields_item_private.h"

#include "psc_bnd.h"
#include "psc_fields_as_c.h"

#include "libpsc/psc_moments/common_moments.c"

static void
do_1st_calc_kh(int p, fields_t *pf, particles_t *pp)
{
  particle_real_t fnqs = sqr(ppsc->coeff.alpha) * ppsc->coeff.cori / ppsc->coeff.eta;
  particle_real_t dxi = 1.f / ppsc->dx[0], dyi = 1.f / ppsc->dx[1], dzi = 1.f / ppsc->dx[2];

  struct psc_patch *patch = &ppsc->patch[p];
  for (int n = 0; n < pp->n_part; n++) {
    particle_t *part = particles_get_one(pp, n);
    int m = particle_single_kind(part);
    DEPOSIT_TO_GRID_1ST_CC(part, pf, m, particle_qni(part));
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
  psc_mfields_zero_range(res, 0, 4);
  
  psc_foreach_patch(ppsc, p) {
    do_1st_calc_kh(p, psc_mfields_get_patch_c(res, p),
		   psc_mparticles_get_patch(particles, p));
  }
  prof_stop(pr);

  psc_mparticles_put_cf(particles, particles_base, MP_DONT_COPY);

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
  double pert;
  double pert_vpic;

  // calculated from the above
  double B0;
  double v0z;
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
  { "pert"          , VAR(pert)            , PARAM_DOUBLE(.0)            },
  { "pert_vpic"     , VAR(pert_vpic)       , PARAM_DOUBLE(.0)            },
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

  struct psc_kind kinds[NR_KH_KINDS] = {
    [KH_ELECTRON1] = { .q = -1., .m = 1, },
    [KH_ELECTRON2] = { .q = -1., .m = 1, },
    [KH_ION1]      = { .q =  1.,         },
    [KH_ION2]      = { .q =  1.,         },
  };
  psc_set_kinds(psc, NR_KH_KINDS, kinds);

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

  psc->domain.bnd_part_lo[0] = BND_PART_PERIODIC;
  psc->domain.bnd_part_hi[0] = BND_PART_PERIODIC;
  psc->domain.bnd_part_lo[1] = BND_PART_REFLECTING;
  psc->domain.bnd_part_hi[1] = BND_PART_REFLECTING;
  psc->domain.bnd_part_lo[2] = BND_PART_PERIODIC;
  psc->domain.bnd_part_hi[2] = BND_PART_PERIODIC;

  // FIXME: can only use 1st order pushers with current conducting wall b.c.
  psc_push_particles_set_type(psc->push_particles, "1vb");
  psc_moments_set_type(psc->moments, "c_1st_cc");
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
  double mi = me * kh->mi_over_me;
  double B0 = sqrt(me) / (kh->wpe_over_wce);
  double vAe = B0 / sqrt(me);
  //  double vAi = B0 / sqrt(mi);
  double vAe_plane = vAe * cos(kh->theta_V);
  //  double v0 = .5 * vAi;
  double v0 = 2 * vAe_plane;
  double v0z = v0 * cos(kh->theta_V - kh->theta_B);
  double Te = kh->beta * (1. / (1. + kh->Ti_over_Te)) * sqr(B0) / 2.;
  double Ti = kh->beta * (1. / (1. + 1./kh->Ti_over_Te)) * sqr(B0) / 2.;
  mpi_printf(MPI_COMM_WORLD, "psc/kh: v0=%g v0z=%g\nvAe=%g vAe_plane=%g\n",
	     v0, v0z, vAe, vAe_plane);
  mpi_printf(MPI_COMM_WORLD, "psc/kh: Te %g Ti %g\n", Te, Ti);
  mpi_printf(MPI_COMM_WORLD, "psc/kh: lambda_De %g\n", sqrt(Te));

  kh->B0 = B0;
  kh->v0z = v0z;

  // set particle kind parameters
  assert(psc->prm.nr_kinds == NR_KH_KINDS);
  psc->kinds[KH_ELECTRON1].T = Te;
  psc->kinds[KH_ELECTRON2].T = Te;
  psc->kinds[KH_ION1].m = mi;
  psc->kinds[KH_ION2].m = mi;
  psc->kinds[KH_ION1].T = Ti;
  psc->kinds[KH_ION2].T = Ti;

  psc_setup_super(psc);
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

  double yl = psc->domain.length[1], zl = psc->domain.length[2];
  double vz = kh->v0z * tanh((x[1] - .5 * yl * (1. + kh->pert * sin(2*M_PI * x[2] / zl))) / kh->delta);
  vz += kh->pert_vpic * kh->v0z * sin(.5 * x[2] / kh->delta) * exp(-sqr(x[1] - .5 * yl)/sqr(kh->delta));

  npt->p[2] = vz;
  switch (kind) {
  case KH_ELECTRON1:
  case KH_ION1:
    *wni = 1.;
    if (vz < 0.) {
      npt->n = 0.;
    } else {
      npt->n = 1.;
    }
    break;
  case KH_ELECTRON2:
  case KH_ION2:
    *wni = 1. + 1e-6;
    if (vz < 0.) {
      npt->n = 1.;
    } else {
      npt->n = 0.;
    }
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

  mparticles_c_t *particles = NULL;
  if (!count_only) {
    particles = psc_mparticles_get_c(psc->particles, MP_DONT_COPY);
  }

  psc_foreach_patch(psc, p) {
    particles_c_t *pp = NULL;
    if (!count_only) {
      pp = psc_mparticles_get_patch_c(particles, p);
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

	    struct psc_particle_npt npt = {
	      .q    = psc->kinds[kind].q,
	      .m    = psc->kinds[kind].m,
	      .n    = psc->kinds[kind].n,
	      .T[0] = psc->kinds[kind].T,
	      .T[1] = psc->kinds[kind].T,
	      .T[2] = psc->kinds[kind].T,
	      // rest is initialized to zero
	    };
	    double wni;
	    psc_kh_init_npt(psc, kind, xx, &npt, &wni);
	    
	    int n_in_cell = npt.n * psc->prm.nicell + .5;
	    if (count_only) {
	      i += n_in_cell;
	    }
	    if (!count_only) {
	      for (int cnt = 0; cnt < n_in_cell; cnt++) {
		particle_c_t *p = particles_c_get_one(pp, i++);
		
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
    psc_mparticles_put_c(particles, psc->particles, 0);
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
