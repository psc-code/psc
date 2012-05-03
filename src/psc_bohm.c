
#include <psc.h>
#include <psc_push_particles.h>
#include <psc_moments.h>
#include <psc_push_fields.h>
#include <psc_sort.h>
#include <psc_balance.h>
#include <psc_particles_as_c.h>
#include <psc_event_generator_private.h>

#include <mrc_params.h>
#include <mrc_profile.h>

#include <math.h>
#include <time.h>

struct psc_es1 {
  // parameters
  double vth_e;
  double vth_i;
  double mi_over_me;
};

#define to_psc_es1(psc) mrc_to_subobj(psc, struct psc_es1)

#define VAR(x) (void *)offsetof(struct psc_es1, x)
static struct param psc_es1_descr[] = {
  { "vth_e"         , VAR(vth_e)            , PARAM_DOUBLE(0.2)           },
  { "vth_i"         , VAR(vth_i)            , PARAM_DOUBLE(0.01)          },
  { "mi_over_me"    , VAR(mi_over_me)       , PARAM_DOUBLE(100)           },

  {},
};
#undef VAR

// ----------------------------------------------------------------------
// psc_es1_create

static void
psc_es1_create(struct psc *psc)
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

  psc->prm.nr_kinds = 2;
  psc->prm.nicell = 50;
  psc->prm.cfl = 0.98;

  psc->domain.length[0] = 1.; // no x-dependence
  psc->domain.length[1] = 1.; // no y-dependence
  psc->domain.length[2] = 20.;

  psc->domain.gdims[0] = 1;
  psc->domain.gdims[1] = 1;
  psc->domain.gdims[2] = 100;

  psc->domain.bnd_fld_lo[0] = BND_FLD_PERIODIC;
  psc->domain.bnd_fld_hi[0] = BND_FLD_PERIODIC;
  psc->domain.bnd_fld_lo[1] = BND_FLD_PERIODIC;
  psc->domain.bnd_fld_hi[1] = BND_FLD_PERIODIC;
  psc->domain.bnd_fld_lo[2] = BND_FLD_CONDUCTING_WALL;
  psc->domain.bnd_fld_hi[2] = BND_FLD_CONDUCTING_WALL;
  psc->domain.bnd_part_lo[0] = BND_PART_PERIODIC;
  psc->domain.bnd_part_hi[0] = BND_PART_PERIODIC;
  psc->domain.bnd_part_lo[1] = BND_PART_PERIODIC;
  psc->domain.bnd_part_hi[1] = BND_PART_PERIODIC;
  psc->domain.bnd_part_lo[2] = BND_PART_REFLECTING;
  psc->domain.bnd_part_hi[2] = BND_PART_ABSORBING;

  psc_moments_set_type(psc->moments, "1st_cc");
}

// ----------------------------------------------------------------------
// psc_es1_init_field

static double
psc_es1_init_field(struct psc *psc, double x[3], int m)
{
  //struct psc_es1 *es1 = to_psc_es1(psc);


  switch (m) {
 
  default: return 0.;
  }
}

// ----------------------------------------------------------------------
// psc_es1_setup_particles

static void
psc_es1_init_npt(struct psc *psc, int kind, double x[3],
		struct psc_particle_npt *npt, double *wni)
{
  struct psc_es1 *es1 = to_psc_es1(psc);

  npt->n = 1.;

  switch (kind) {
  case 0: // electrons
    npt->q = -1.;
    npt->m = 1.;
    npt->T[0] = sqr(es1->vth_e);
    npt->T[1] = sqr(es1->vth_e);
    npt->T[2] = sqr(es1->vth_e);
    break;
  case 1: // ions
    npt->q = 1.;
    npt->m = es1->mi_over_me;
    npt->T[0] = sqr(es1->vth_i)*npt->m;
    npt->T[1] = sqr(es1->vth_i)*npt->m;
    npt->T[2] = sqr(es1->vth_i)*npt->m;
    break;
  default:
    assert(0);
  }
}


// ======================================================================
// psc_es1_ops

struct psc_ops psc_es1_ops = {
  .name             = "es1",
  .size             = sizeof(struct psc_es1),
  .param_descr      = psc_es1_descr,
  .create           = psc_es1_create,
  .init_field       = psc_es1_init_field,
  .init_npt         = psc_es1_init_npt,
};

// ======================================================================
// particle seeding

static void
seed_patch(struct psc *psc, int p, particles_t *pp)
{
  //  struct psc_es1 *es1 = to_psc_es1(psc);

  float r = random() / (float) RAND_MAX;

  particles_realloc(pp, pp->n_part + 2);

  double xi = CRDX(p, 0);
  double yi = CRDY(p, 0);
  double zi = r * psc->domain.length[2];

  struct psc_particle_npt npt = {};
  particle_t *prt;

  // electrons
  prt = particles_get_one(pp, pp->n_part++);
  prt->xi = xi;
  prt->yi = yi;
  prt->zi = zi;
  prt->wni = 1.;
  npt.q = -1.;
  npt.m = 1.;
  npt.T[0] = 0.01;
  npt.T[1] = 0.01;
  npt.T[2] = 0.01;
  psc_setup_particle(psc, prt, &npt);

  // ions
  prt = particles_get_one(pp, pp->n_part++);
  prt->xi = xi;
  prt->yi = yi;
  prt->zi = zi;
  prt->wni = 1.;
  npt.q = 1.;
  npt.m = 100.;
  npt.T[0] = 0.01;
  npt.T[1] = 0.01;
  npt.T[2] = 0.01;
  psc_setup_particle(psc, prt, &npt);
}

void
psc_event_generator_bohm_run(struct psc_event_generator *gen,
			     mparticles_base_t *mparticles, mfields_base_t *mflds,
			     mphotons_t *mphotons)
{
  psc_foreach_patch(ppsc, p) {
    seed_patch(ppsc, p, psc_mparticles_get_patch(mparticles, p));
  }
}

// ======================================================================
// psc_event_generator: subclass "bohm"

struct psc_event_generator_ops psc_event_generator_bohm_ops = {
  .name                  = "bohm",
  .run                   = psc_event_generator_bohm_run,
};

// ======================================================================
// main

int
main(int argc, char **argv)
{
  mrc_class_register_subclass(&mrc_class_psc_event_generator,
			      &psc_event_generator_bohm_ops);

  return psc_main(&argc, &argv, &psc_es1_ops);
}
