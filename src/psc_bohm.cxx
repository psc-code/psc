
#include <psc.h>
#include <psc_push_particles.h>
#include <psc_push_fields.h>
#include <psc_sort.h>
#include <psc_balance.h>
#include <psc_particles_as_double.h>
#include <psc_event_generator_private.h>

#include <mrc_params.h>
#include <mrc_profile.h>

#include <math.h>
#include <time.h>

struct psc_bohm {
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

#define to_psc_bohm(psc) mrc_to_subobj(psc, struct psc_bohm)

#define VAR(x) (void *)offsetof(struct psc_bohm, x)
static struct param psc_bohm_descr[] = {
  { "Te"            , VAR(Te)            , PARAM_DOUBLE(1)           },
  { "Ti"            , VAR(Ti)            , PARAM_DOUBLE(0.023)       },
  { "mi_over_me"    , VAR(mi_over_me)    , PARAM_DOUBLE(100)         },
  { "L"             , VAR(L)             , PARAM_DOUBLE(0.01)        },
  { "L_source"      , VAR(L_source)      , PARAM_DOUBLE(0.01)        },

  {},
};
#undef VAR


// ----------------------------------------------------------------------
// psc_bohm_create

static void
psc_bohm_create(struct psc *psc)
{
  psc_default_dimensionless(psc);

  psc->prm.nmax = 16000;
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
// psc_bohm_setup

static void
psc_bohm_setup(struct psc *psc)
{
  struct psc_bohm *bohm = to_psc_bohm(psc);
  double M=9.11e-31;
  double C=6.0e7;
  double e=1.6e-19;
  double eps_o=8.85e-12;
  double no=1e15;
  
  bohm->Te_=e*bohm->Te/(M*C*C);
  bohm->Ti_=e*bohm->Ti/(M*C*C);
  double vte=sqrt(2*e*bohm->Te/M);
  double lde=sqrt(eps_o*e*bohm->Te/(no*e*e));
  double de=C*sqrt(eps_o*M/(e*e*no));
  double vti=sqrt(2*e*bohm->Ti/(M*bohm->mi_over_me));
  double cs=sqrt(e*bohm->Te/(M*bohm->mi_over_me));
 
  bohm->L_source_=bohm->L_source/de;
  bohm->cs_=cs/C;
  
  psc->domain.length[2] = bohm->L/de;

  psc_setup_super(psc);

  printf("lambda_de=%g(%g) dz=%g\n", sqrt(bohm->Te_), lde, psc->patch[0].dx[2]);
  printf("v_te=%g(%g)\n", vte/C, vte);
  printf("v_ti=%g(%g)\n", vti/C, vti);
  printf("cs=%g(%g)\n", cs/C, cs);
  printf("de=%g\n", de);
  printf("L=%g(%g)\n", bohm->L/de, bohm->L);

}

// ----------------------------------------------------------------------
// psc_bohm_init_field

static double
psc_bohm_init_field(struct psc *psc, double x[3], int m)
{
  //struct psc_bohm *bohm = to_psc_bohm(psc);


  switch (m) {
 
  default: return 0.;
  }
}

// ----------------------------------------------------------------------
// psc_bohm_init_npt

static void
psc_bohm_init_npt(struct psc *psc, int kind, double x[3],
		struct psc_particle_npt *npt)
{
  struct psc_bohm *bohm = to_psc_bohm(psc);

  npt->n = 1.;

  switch (kind) {
  case 0: // electrons
    npt->q = -1.;
    npt->m = 1.;
    npt->T[0] = bohm->Te_;
    npt->T[1] = bohm->Te_;
    npt->T[2] = bohm->Te_;
    break;
  case 1: // ions
    npt->q = 1.;
    npt->m = bohm->mi_over_me;
    npt->T[0] = bohm->Ti_;
    npt->T[1] = bohm->Ti_;
    npt->T[2] = bohm->Ti_;
    break;
  default:
    assert(0);
  }
}


// ======================================================================
// psc_bohm_ops

struct psc_ops psc_bohm_ops = {
  .name             = "bohm",
  .size             = sizeof(struct psc_bohm),
  .param_descr      = psc_bohm_descr,
  .create           = psc_bohm_create,
  .setup            = psc_bohm_setup,
  .init_field       = psc_bohm_init_field,
  .init_npt         = psc_bohm_init_npt,
};

// ======================================================================
// particle seeding

static void
seed_patch(struct psc *psc, struct psc_mparticles *mprts, int p)
{
  struct psc_bohm *bohm = to_psc_bohm(psc);
  //particle_range_t prts = mparticles_t(mprts)[p].range();
  
  psc_foreach_3d(psc, p, ix, iy, iz, 0, 0) {

    double Iono_rate;
    double xx[3] = { CRDX(p, ix), CRDY(p, iy), .5 * (CRDZ(p, iz) + CRDZ(p, iz+1)) };

    if (xx[2] < bohm->L_source_) {
    // Number of new particles created per unit time and cell
      Iono_rate = psc->prm.nicell * 0.6 * bohm->cs_ / bohm->L_source_;
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

    mparticles_t(mprts)[p].reserve(mparticles_t(mprts)[p].size() + 2*N_new);

    assert(0);
#if 0
    struct psc_particle_npt npt = {};

    for (int n = 0; n < N_new; n++) {
      particle_t prt;

      // electrons
      prt.wni = 1.;
      npt.q = -1.;
      npt.m = 1.;
      npt.T[0] = bohm->Te_;
      npt.T[1] = bohm->Te_;
      npt.T[2] = bohm->Te_;
      npt.kind = KIND_ELECTRON;
      psc_setup_particle(psc, &prt, &npt, p, xx);
      mprts[p].push_back(prt);

      // ions
      prt.wni = 1.;
      npt.q = 1.;
      npt.m = bohm->mi_over_me;
      npt.T[0] = bohm->Ti_;
      npt.T[1] = bohm->Ti_;
      npt.T[2] = bohm->Ti_;
      npt.kind = KIND_ION;
      psc_setup_particle(psc, &prt, &npt, p, xx);
      mprts[p].push_back(prt);
    }
#endif
  } foreach_3d_end;
}

void
psc_event_generator_bohm_run(struct psc_event_generator *gen,
			     struct psc_mparticles *mprts, struct psc_mfields *mflds)
{
  psc_foreach_patch(ppsc, p) {
    seed_patch(ppsc, mprts, p);
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

  return psc_main(&argc, &argv, &psc_bohm_ops);
}
