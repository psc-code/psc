
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
  double Te;
  double Ti;
  double mi_over_me;
  double L;
  double L_source;

  double Te_;
  double Ti_;
  double cs_;
  double L_source_;
};

#define to_psc_es1(psc) mrc_to_subobj(psc, struct psc_es1)

#define VAR(x) (void *)offsetof(struct psc_es1, x)
static struct param psc_es1_descr[] = {
  { "Te"            , VAR(Te)            , PARAM_DOUBLE(1)           },
  { "Ti"            , VAR(Ti)            , PARAM_DOUBLE(0.023)       },
  { "mi_over_me"    , VAR(mi_over_me)    , PARAM_DOUBLE(100)         },
  { "L"             , VAR(L)             , PARAM_DOUBLE(0.01)        },
  { "L_source"      , VAR(L_source)      , PARAM_DOUBLE(0.01)        },

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

  // psc_moments_set_type(psc->moments, "c_1st_cc");
}

// ----------------------------------------------------------------------
// psc_es1_setup

static void
psc_es1_setup(struct psc *psc)
{
  struct psc_es1 *es1 = to_psc_es1(psc);
  double M=9.11e-31;
  double C=6.0e7;
  double e=1.6e-19;
  double eps_o=8.85e-12;
  double no=1e15;
  
  es1->Te_=e*es1->Te/(M*C*C);
  es1->Ti_=e*es1->Ti/(M*C*C);
  double vte=sqrt(2*e*es1->Te/M);
  double lde=sqrt(eps_o*e*es1->Te/(no*e*e));
  double de=C*sqrt(eps_o*M/(e*e*no));
  double vti=sqrt(2*e*es1->Ti/(M*es1->mi_over_me));
  double cs=sqrt(e*es1->Te/(M*es1->mi_over_me));
 
  es1->L_source_=es1->L_source/de;
  es1->cs_=cs/C;
  
  psc->domain.length[2] = es1->L/de;

  psc_setup_super(psc);

  printf("lambda_de=%g(%g) dz=%g\n", sqrt(es1->Te_), lde, psc->dx[2]);
  printf("v_te=%g(%g)\n", vte/C, vte);
  printf("v_ti=%g(%g)\n", vti/C, vti);
  printf("cs=%g(%g)\n", cs/C, cs);
  printf("de=%g\n", de);
  printf("L=%g(%g)\n", es1->L/de, es1->L);

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
// psc_es1_init_npt

static void
psc_es1_init_npt(struct psc *psc, int kind, double x[3],
		struct psc_particle_npt *npt)
{
  struct psc_es1 *es1 = to_psc_es1(psc);

  npt->n = 1.;

  switch (kind) {
  case 0: // electrons
    npt->q = -1.;
    npt->m = 1.;
    npt->T[0] = es1->Te_;
    npt->T[1] = es1->Te_;
    npt->T[2] = es1->Te_;
    break;
  case 1: // ions
    npt->q = 1.;
    npt->m = es1->mi_over_me;
    npt->T[0] = es1->Ti_;
    npt->T[1] = es1->Ti_;
    npt->T[2] = es1->Ti_;
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
  .setup            = psc_es1_setup,
  .init_field       = psc_es1_init_field,
  .init_npt         = psc_es1_init_npt,
};

// ======================================================================
// particle seeding

static void
seed_patch(struct psc *psc, int p, particles_t *pp)
{
  struct psc_es1 *es1 = to_psc_es1(psc);
  
  psc_foreach_3d(psc, p, ix, iy, iz, 0, 0) {

    double Iono_rate;
    double x = CRDX(p, ix);
    double y = CRDY(p, iy);
    double z = CRDZ(p, iz);

    if (z < es1->L_source_) {
    // Number of new particles created per unit time and cell
      Iono_rate = psc->prm.nicell * 0.6 * es1->cs_ / es1->L_source_;
    } else {
      Iono_rate = 0;
    }

    double Iono_count = Iono_rate*psc->dt;


    // inverse Poisson, hope Sdt is not huge!

    float r = random() / (float) RAND_MAX;
    int N_new = 0;
    float f = 1;
    float q = r * exp(Iono_count) - 1;
    while (q > 0) {
      N_new++;
      f*= Iono_count/N_new;
      q-= f;
    }

    particles_realloc(pp, pp->n_part + 2*N_new);

    struct psc_particle_npt npt = {};
    particle_t *prt;

    for (int n=0; n < N_new; n++) {

      // electrons
      prt = particles_get_one(pp, pp->n_part++);
      prt->xi = x;
      prt->yi = y;
      prt->zi = z;
      prt->wni = 1.;
      npt.q = -1.;
      npt.m = 1.;
      npt.T[0] = es1->Te_;
      npt.T[1] = es1->Te_;
      npt.T[2] = es1->Te_;
      psc_setup_particle(psc, prt, 0, &npt);
      
      // ions
      prt = particles_get_one(pp, pp->n_part++);
      prt->xi = x;
      prt->yi = y;
      prt->zi = z;
      prt->wni = 1.;
      npt.q = 1.;
      npt.m = es1->mi_over_me;
      npt.T[0] = es1->Ti_;
      npt.T[1] = es1->Ti_;
      npt.T[2] = es1->Ti_;
      psc_setup_particle(psc, prt, 1, &npt);
      
    }
  } foreach_3d_end;

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
